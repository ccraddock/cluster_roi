#!/usr/bin/env python
import argparse
import os
import subprocess
import json
import sys

import datetime
import time
import urllib.request

import pyclusterroi as pc


def run(command, env={}):
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                               shell=True, env=env)
    while True:
        line = process.stdout.readline()
        line = str(line)[:-1]
        print(line)
        if line == '' and process.poll() is not None:
            break


def pyclusterroi(input_arguments):
    """
    functionality for command line interface, we follow the bids app spec so that it can be used there too

    :param input_arguments: command line argument
    :return: 0 on success, 1 on error
    """

    parser = argparse.ArgumentParser(description='PyClusterROI {0} Command Line Tool'.format(pc.__version__))

    parser.add_argument('bids_dir',
                        help='This directory should have a /derivatives subdirectory were the input 4D preprocessed '
                             'fMRI data or results from a participant-level clustering have been organized according to '
                             'the BIDS derivatives standard. be either 4D preprocessed fMRI data or results from a '
                             'first-level clustering that should be located in the /derivatives subdirectory of the '
                             'bids directory. Use the format s3://bucket/path/to/bidsdir to read data directly from an '
                             'S3 bucket. This may require AWS S3 credentials specified via the --aws_input_credentials '
                             'option.')

    parser.add_argument('output_dir',
                        help='The directory where the output files should be stored. Output files will be arranged '
                             'according to the BIDS derivative standard in a /derivative subdirectory. Use the format '
                             's3://bucket/path/to/bidsdir to write data directly to an S3 bucket. This may (will '
                             'likely) require AWS S3 credentials specified via the --aws_output_credentials option.')

    parser.add_argument('analysis_level',
                        help='Level of the analysis that will  be performed. Multiple participant level analyses can '
                             'be run independently (in parallel) using the same output_dir. Group level analyses can '
                             'be accomplished by performing a second-level parcellation of results from participant '
                             'specific parcellations (see --2level flag) or by parcellating the average similarity '
                             'matrix calculated from the functional data (see --gmean flag). test_participant will '
                             'and test_group run through the entire configuration process for the respective analysis '
                             'but will not execute the pipeline.',
                        choices=['participant', 'group', 'test_participant', 'test_group'])

    parser.add_argument('--similarity_metric',
                        help='The method used to calculate similarity between neighboring voxels. Possible values '
                             'are "tcorr": the Pearson correlation between voxel time series, "scorr": Pearson '
                             'correlation between whole brain functional connectivity maps generated from the '
                             'respective voxel time series, "ones": 1 if voxels are neighbors, 0 otherwise, '
                             '"cluster": 1 if voxels are in the same cluster region and 0 if not. "cluster" is '
                             'intended for second level clustering where the inputs are individual level parcellation '
                             'results. Defaults to "tcorr"',
                        choices=['tcorr', 'scorr', 'ones', 'cluster'])

    parser.add_argument('--group_level_method',
                        help='The method used to perform group level clustering. This can either be "gmean" in which a '
                             'group similarity matrix is constructed by averaging the similarity matrices calculated '
                             'from each group member\'s 4D date and then submitted to clustering. For the "2level" '
                             'method parcellation results for each participant are transformed into a similarity '
                             'matrix by setting the similarity between voxels in the same parcellation region to "1" '
                             'and "0" elsewhere, these similarity matrices are avaeraged across participants and then '
                             'submitted to clustering. For "mean" the derivative type will default to "preproc", for '
                             '"2level" it will default to "cluster". Defaults to "2level"',
                        choices=['gmean', '2level'])

    parser.add_argument('--group_name',
                        help='Consistent with the BIDS format for participant labelling, group results will be '
                             'labelled with grp-<group_name> and a json file will be generated listing all of '
                             'the datasets that went into the group. If not given, a date time stamp will be used.',
                        default=None)

    parser.add_argument('--number_of_clusters',
                        help='The number of regions desired from the parcellation. Several different levels can be '
                             'specified in a comma seperated list with no intervening spaces (e.g. 10,20,30,40,50). '
                             'Can also specify with using start:end:stride notation (e.g. 10:60:10 is equivalent to '
                             '10,20,30,40,50). This is a required parameter, and must be specified either on the '
                             'command line or in the configuration in file.')

    parser.add_argument('--thresh',
                        help='When calculating the similarity between voxels, values below this threshold are set to '
                             '0. This is intended to remove values that may have arisen due to chance. A suitable '
                             'choice can be calculated from a phase randomization monte carlo simulation of two '
                             'randomly chosen voxels. Default = 0.5.',
                        default='0.5')

    parser.add_argument('--mask',
                        help='Path to mask to be used to restrict the parcellation - typically a grey matter mask. '
                             'This is a required parameter, and must be specified either on the command line or in '
                             'the configuration in file.')

    parser.add_argument('--config_file',
                        help='A .json file containing configuration parameters. Provides more flexibility for '
                             'determining which data in the bids directory is processed and other parameters than is '
                             'available with the command line. A configuration file is written into the output '
                             'directory at every execution and can be customized for future runs. The "test_config" '
                             'analysis_level will produce a configuration file based on the command line flags without '
                             'executing the proscribed parcellation. Command line parameters take precedence over '
                             'parameters in the config file.')

    parser.add_argument('--data_config_file',
                        help='json file containing a list of paths to the data to be used in the parcellation. '
                             'Eliminates the overhead of re-reading the bids input directory at every execution, which '
                             'is useful for batch processing (see --participant_ndx).')

    parser.add_argument('--aws_input_credentials',
                        help='Credentials for reading from S3. If not provided and s3 paths are specified in the data '
                             'config we will try to access the bucket anonymously',
                        default=None)

    parser.add_argument('--aws_output_credentials',
                        help='Credentials for writing to S3. If not provided and s3 paths are specified in the output '
                             'directory we will try to access the bucket anonymously',
                        default=None)

    parser.add_argument('--n_cpus',
                        help='Number of threads to use for the parcellation. This specifically parallelizes the '
                             'calculation of each voxel\'s connectivity, and the binarization of clustering results '
                             'across different clustering levels. For smaller datasets (4mm voxels), this really only '
                             'improves the performance for scorr connectivity based analyses. Will likely be more '
                             'significant at finer resolution data.',
                        default='1')

    parser.add_argument('--participant_label',
                        help='The label of the participant that should be analyzed. The label corresponds to '
                             'sub-<participant_label> from the BIDS spec (so it does not include "sub-"). If this '
                             'parameter is not provided all participants should be analyzed. Multiple participants can '
                             'be specified with a space separated list. To work correctly this should come at the end '
                             'of the command line.',
                        nargs="+")

    parser.add_argument('--participant_ndx',
                        help='The index of the participant that should be analyzed. This corresponds to the index of '
                             'the participant in the data list file. This was added to make it easier to accommodate '
                             'cluster array jobs. Only a single participant will be analyzed. Can be used with '
                             'participant label, in which case it is the index into the list that follows the '
                             'participant_label flag.',
                        default=None)

    parser.add_argument('--pipeline_name',
                        help='The name of the pipeline that generated the files to be parcellated. Only needed if '
                             'there is more than one pipeline subdirectory in the bids derivatives file structure '
                             'containing derivatives that match the other input criteria (task, derivative, etc) and '
                             'you only want to analyze data from one of them.',
                        default=None)

    parser.add_argument('--session',
                        help='The BIDS session value for the data to be parcellated. Defaults to all sessions.',
                        default=None)

    parser.add_argument('--consolidate_sessions',
                        help='Consolidate data from different sessions into single parcellation by '
                             'averaging connectivity matrices calculated for each session.',
                        action='store_true')

    parser.add_argument('--task',
                        help='The BIDS task value for the data to be parcellated. Defaults to "rest".',
                        default='rest')

    parser.add_argument('--consolidate_tasks',
                        help='Consolidate data from different tasks into single parcellation by averaging connectivity '
                             'matrices calculated for each task.',
                        action='store_true')

    parser.add_argument('--run',
                        help='The run index of the data to be parcellated. Only need if there is more than one run '
                             'that matches the other input criteria (task, derivative, etc) and you only want to '
                             'parcellate one of them. Default behavior is to seperately parcellate data from each run. '
                             'use --consolidate_runs to change this behavior.',
                        default=None)

    parser.add_argument('--consolidate_runs',
                        help='Consolidate data from different runs into single parcellation by averaging connectivity '
                             'matrices calculated for each run.',
                        action='store_true')

    parser.add_argument('--derivative',
                        help='The BIDS derivative type for the data to be parcellated. Defaults to "preproc".',
                        default='preproc')

    parser.add_argument('--variant',
                        help='The BIDS variant type for the data to be parcellated. Only need if there is more than '
                             'variant that matches the other input criteria (task, derivative, etc) and you only want '
                             'to parcellate one of them.',
                        default=None)

    parser.add_argument('--bids_validator_config',
                        help='JSON file specifying configuration of bids-validator: See '
                             'https://github.com/INCF/bids-validator for more info')

    parser.add_argument('--skip_bids_validator',
                        help='Skips BIDS validation.',
                        action='store_true')

    parser.add_argument('--debug',
                        help='Print out additional information for debugging.',
                        action='store_true')

    # get the command line arguments
    args = parser.parse_args(input_arguments)

    # create a timestamp for writing config files
    timestamp = datetime.datetime.fromtimestamp(time.time()).strftime('%Y%m%d%H%M%S')

    debug_flag = False

    if args.debug is True:
        debug_flag = True

    if debug_flag is True:
        print(input_arguments)
        print(args)

    # check to make sure that the input directory exists
    if not args.bids_dir.lower().startswith("s3://") and not os.path.exists(args.bids_dir):
        print("Error! Could not find {0}".format(args.bids_dir))
        sys.exit(0)

    # check to make sure that the output directory exists
    if not args.output_dir.lower().startswith("s3://") and not os.path.exists(args.output_dir):
        print("Error! Could not find {0}".format(args.output_dir))
        sys.exit(0)

    parameters_dict = {'pipeline': None,
                       'derivative': ['preproc'],
                       'similarity_metric': 'tcorr',
                       'group_level_method': '2level',
                       'session': None,
                       'consolidate_sessions': False,
                       'task': ['rest'],
                       'consolidate_tasks': False,
                       'run': None,
                       'consolidate_runs': False,
                       'variant': None,
                       'number_of_threads': 1,
                       'number_of_clusters': [],
                       'mask': None,
                       'aws_input_credentials': None,
                       'aws_output_credentials': None,
                       'participant_label': None,
                       'participant_index': None}

    # otherwise, if we are running group, participant, or dry run we
    # begin by conforming the configuration
    if args.config_file:
        if args.config_file.lower().startswith("s3://"):
            with urllib.request.urlopen('http://python.org/') as response:
                parameters_dict.update(json.loads(response.read()))
        else:
            parameters_dict.update(json.load(open(os.path.realpath(args.config_file), 'r')))

    if args.pipeline_name:
        parameters_dict['pipeline_name'] = args.pipeline.split(',')

    if args.session:
        parameters_dict['session'] = args.session.split(',')

    if args.consolidate_sessions:
        parameters_dict['consolidate_sessions'] = args.consolidate_sessions

    if args.derivative:
        parameters_dict['derivative'] = args.derivative.split(',')

    if args.task:
        parameters_dict['task'] = args.task.split(',')

    if args.consolidate_tasks:
        parameters_dict['consolidate_tasks'] = args.consolidate_tasks

    if args.variant:
        parameters_dict['variant'] = args.variant.split(',')

    if args.run:
        parameters_dict['run'] = args.run.split(',')

    if args.consolidate_runs:
        parameters_dict['consolidate_runs'] = args.consolidate_runs

    if args.similarity_metric:

        if 'cluster' in args.similarity_metric:
            if 'participant' in args.analysis_level:
                raise ValueError('Participant level analysis cannot be performed with the cluster similarity metric.')
            if 'roi' not in parameters_dict['derivative']:
                print('The cluster similarity only works with the "roi" derivative type data. Making the changing the '
                      'derivative type to "roi" ...')
                parameters_dict['derivative'] = 'roi'
            if args.group_level_method is 'gmean':
                print('Cannot use "gmean" group level method with "cluster" similarity. Changing to 2level')
                parameters_dict['group_level_method'] = '2level'

        parameters_dict['similarity_metric'] = args.similarity_metric

    if args.group_level_method:
        parameters_dict['group_level_method'] = args.group_level_method

    if args.analysis_level is 'group':
        
        if parameters_dict['group_level_method'] is '2level':

            if parameters_dict['similarity_metric'] is not 'cluster':
                print('The 2corr method for group level analysis requires the cluster similarity metric. Changing ...')
                parameters_dict['similarity_metric'] = 'cluster'
            if parameters_dict['derivative'] is not 'roi':
                print('The 2corr method requires derivative type to be "roi". Changing ... ')
                parameters_dict['derivative'] = 'roi'

        elif args.group_level_method is 'gmean':
            if parameters_dict['similarity_metric'] is 'cluster':
                raise ValueError('The gmean method for group level analysis cannot use the cluster similarity metric.')
            if parameters_dict['derivative'] is 'roi':
                raise ValueError('The gmean method for group level analysis cannot use derivative type "roi".')

        if args.group_name:
            parameters_dict['group_name'] = args.group_name
        else:
            parameters_dict['group_name'] = timestamp

    if args.number_of_clusters:
        if ',' in args.number_of_clusters:
            parameters_dict["number_of_clusters"] = []
            for k in args.number_of_clusters.split(','):
                parameters_dict["number_of_clusters"].append(int(k))

        elif ':' in args.number_of_clusters:
            string_values = args.number_of_clusters.split(':')
            if len(string_values) == 3:
                parameters_dict["number_of_clusters"] = list(range(int(string_values[0]),
                                                                   int(string_values[1]),
                                                                   int(string_values[2])))
            else:
                raise ValueError("Could not process {0}, is it correctly formatted?".format(args.number_of_clusters))

        else:
            parameters_dict["number_of_clusters"] = int(args.number_of_clusters)

    if args.thresh:
        parameters_dict['thresh'] = float(args.thresh)
        if 0 > parameters_dict['thresh'] or parameters_dict['thresh'] > 1:
            raise ValueError(
                'Inappropriate choice for threshhold ({0}), should be [0,1)'.format(parameters_dict['thresh']))

    if args.mask:
        if not os.path.isfile(args.mask):
            raise ValueError("Could not find mask file {0}.".format(args.mask))
        parameters_dict["mask"] = args.mask
    else:
        raise ValueError("Mask is a required parameter.")

    if args.n_cpus:
        parameters_dict["number_of_threads"] = int(args.n_cpus)

    if args.aws_input_credentials:
        if os.path.isfile(args.aws_input_credentials):
            parameters_dict['aws_input_credentials'] = args.aws_input_credentials
        else:
            raise ValueError("Could not find aws credentials {0}".format(args.aws_input_credentials))

    if args.aws_output_credentials:
        if os.path.isfile(args.aws_output_credentials):
            parameters_dict['aws_output_credentials'] = args.aws_output_credentials
        else:
            raise IOError("Could not find aws credentials {0}".format(args.aws_output_credentials))

    if args.participant_label:
        parameters_dict['participant_label'] = args.participant_label

    if args.participant_ndx:
        parameters_dict['participant_index'] = int(args.participant_ndx)

    print("\n-- Running pyClusterROI with parameters: --")
    print("BIDS input directory: {0}".format(args.bids_dir))
    print("Output directory: {0}".format(args.output_dir))
    print("Configuration file: {0}".format(args.config_file))
    print("Data Configuration file: {0}".format(args.data_config_file))
    print("Mask: {0}".format(parameters_dict['mask']))
    print("Analysis level: {0}".format(args.analysis_level))
    if 'group' in args.analysis_level:
        print("Group level method: {0}".format(parameters_dict['group_level_method']))

    print("Cluster sizes: {0}".format(",".join([str(k) for k in parameters_dict['number_of_clusters']])))
    print("Similarity metric: {0}".format(parameters_dict['similarity_metric']))

    print("Derivative: {0}".format(parameters_dict['derivative']))
    print("Task: {0}".format(parameters_dict['task']))
    print("Run: {0}".format(parameters_dict['run']))
    print("Variant: {0}".format(parameters_dict['variant']))

    print("Number of threads: {0}".format(parameters_dict['number_of_threads']))
    print("AWS input credentials: {0}".format(parameters_dict['aws_input_credentials']))
    print("AWs output credentials: {0}".format(parameters_dict['aws_output_credentials']))

    if parameters_dict['participant_label']:
        print("Participant labels: {0}".format(", ".join(parameters_dict['participant_label'])))

    if parameters_dict['participant_index']:
        print("Participant index: {0}".format(parameters_dict['participant_index']))

    print("--\n")

    # validate input dir (if skip_bids_validator is not set)
    # if args.bids_validator_config:
    #     print("\nRunning BIDS validator with configuration file {config}".format(config=args.bids_validator_config))
    #     run("bids-validator --config {config} {bids_dir}".format(config=args.bids_validator_config,
    #                                                              bids_dir=args.bids_dir))
    # elif args.skip_bids_validator:
    #     print('skipping bids-validator...')
    # else:
    #     print("\nRunning BIDS validator")
    #     run("bids-validator {bids_dir}".format(bids_dir=args.bids_dir))

    # get a list of all of the files matching the criteria defined by the configuration parameters
    bids_root_directory = args.bids_dir.rstrip('//')
    if os.path.basename(bids_root_directory) is not 'derivatives':
        bids_root_directory = os.path.join(args.bids_dir, 'derivatives')

    bids_dictionary = {}
    participant_list = []

    if args.data_config_file:
        with open(args.data_config_file, 'r') as in_stream:
            bids_dictionary = json.load(in_stream)
    # load in data from bids_dir
    else:
        bids_image_filepaths = []
        inclusion_file_count = 0

        if "s3://" not in args.output_dir.lower():
            for path_root, directories, filepaths in os.walk(bids_root_directory):
                for filepath in filepaths:
                    if filepath.endswith(".nii") or filepath.endswith(".nii.gz"):
                        bids_image_filepaths.append(os.path.join(path_root, filepath))
        else:
            raise ValueError("Does not currently support reading data from S3.")

        if debug_flag is True:
            print("Found {0} nifti files from {1}:\n{2}".format(len(bids_image_filepaths), bids_root_directory,
                                                                '\n'.join(bids_image_filepaths)))

        for filepath in bids_image_filepaths:

            filepath_bids_dictionary = {'full_path_name': filepath}

            filepath_parts = filepath.replace(bids_root_directory, '').lstrip('/').split('/')

            if len(filepath_parts) < 3:
                if debug_flag:
                    print("Path has too few components ({0}) to parse, does it conform to bids? Discarding ...".format(
                        filepath_parts))
                continue

            filepath_bids_dictionary['pipeline'] = filepath_parts.pop(0)
            if 'sub-' in filepath_bids_dictionary['pipeline']:
                if debug_flag is True:
                    print("Warning: pipeline name {0} contains sub-, are the directories formatted "
                          "correctly? Discarded".format(filepath_bids_dictionary['pipeline']))
                continue

            if parameters_dict['pipeline'] and \
                    filepath_bids_dictionary['pipeline'] not in parameters_dict['pipeline']:
                if debug_flag is True:
                    print('Discarding {0} because its pipeline {1} is'
                          ' not in the list provided {2}'.format(
                              filepath, filepath_bids_dictionary['pipeline'],
                              ','.join(parameters_dict['pipeline'])))
                continue

            if 'sub-' in filepath_parts[0]:
                filepath_bids_dictionary['participant_label'] = filepath_parts.pop(0).split("-")[1]
                if parameters_dict['participant_label'] and \
                        filepath_bids_dictionary['participant_label'] not in parameters_dict['participant_label']:
                    if debug_flag is True:
                        print('Discarding {0} because its participant label {1} is'
                              ' not in the list provided {2}'.format(
                                  filepath, filepath_bids_dictionary['participant_label'],
                                  ','.join(parameters_dict['participant_label'])))
                    continue
            else:
                if debug_flag is True:
                    print("Did not find participant ID in path {0}, Discarding...".format(filepath_parts[0]))
                continue

            if 'ses-' in filepath_parts[0]:
                filepath_bids_dictionary['ses'] = filepath_parts.pop(0).split("-")[1]
                if parameters_dict['session'] and \
                        filepath_bids_dictionary['session'] not in parameters_dict['session']:
                    if debug_flag is True:
                        print('Discarding {0} because its session {1} is'
                              ' not in the list provided {2}'.format(
                                  filepath, filepath_bids_dictionary['session'],
                                  ','.join(parameters_dict['session'])))
                    continue

            if 'func' not in filepath_parts.pop(0):
                if debug_flag is True:
                    print('Discarding {0} because it does not appear to be a functional file'.format(filepath))
                continue

            filename = filepath_parts.pop(0)
            extension = filename.split(".")[1:]

            if 'nii' not in extension[0]:
                if debug_flag is True:
                    print('Discarding {0} because it does not appear to be a nifti file {1} {2}'.format(filepath,
                                                                                                        filename,
                                                                                                        extension))
                continue

            filename_parts = filename.split(".")[0].split("_")
            for filename_part in filename_parts:
                filename_part_parts = filename_part.split("-")
                if len(filename_part_parts) == 1:
                    if 'filetype' not in filepath_bids_dictionary:
                        filepath_bids_dictionary['filetype'] = filename_part_parts[0]
                    elif 'derivative' not in filepath_bids_dictionary:
                        filepath_bids_dictionary['derivative'] = filename_part_parts[0]
                    else:
                        print('Filename {0} contains more than two singleton segments, '
                              'using last one for derivative'.format(filepath))
                        filepath_bids_dictionary['derivative'] = filename_part_parts[0]
                else:
                    if filename_part_parts[0] in filepath_bids_dictionary:
                        print('Filename {0} appears to have more than one key-value pair for key {1}, using'
                              ' last'.format(filepath, filename_part_parts[0]))
                    else:
                        filepath_bids_dictionary[filename_part_parts[0]] = filename_part_parts[1]

            # verify that the file is appropriate, extreme vetting
            if 'filetype' in filepath_bids_dictionary:
                if 'bold' not in filepath_bids_dictionary['filetype']:
                    if debug_flag is True:
                        print('Discarding {0} because it is not derived from BOLD data'.format(filepath))
                    continue
            else:
                if debug_flag is True:
                    print('Discarding {0} because its filename does not contain a file type'.format(filepath))
                continue

            if 'derivative' in filepath_bids_dictionary:
                if filepath_bids_dictionary['derivative'] not in parameters_dict['derivative']:
                    if debug_flag is True:
                        print('Discarding {0} because its derivative type {1} is not in the list '
                              'provided {2}'.format(filepath, filepath_bids_dictionary['derivative'],
                                                    ','.join(parameters_dict['derivative'])))
                    continue
            else:
                if debug_flag is True:
                    print('Discarding {0} because its filename does not contain a derivative'.format(filepath))
                continue

            if 'task' in filepath_bids_dictionary:
                if filepath_bids_dictionary['task'] not in parameters_dict['task']:
                    if debug_flag is True:
                        print('Discarding {0} because its task {1} is not in the list '
                              'provided {2}'.format(filepath, filepath_bids_dictionary['task'],
                                                    ','.join(parameters_dict['task'])))
                    continue
            else:
                if debug_flag is True:
                    print('Discarding {0} because its filename does not contain task'.format(filepath))
                continue

            if 'variant' in parameters_dict and parameters_dict['variant']:
                if 'variant' in filepath_bids_dictionary:
                    if filepath_bids_dictionary['variant'] not in parameters_dict['variant']:
                        if debug_flag is True:
                            print('Discarding {0} because its variant {1} is not in the list '
                                  'provided {2}'.format(filepath, filepath_bids_dictionary['variant'],
                                                        ','.join(parameters_dict['variant'])))
                        continue
                else:
                    if debug_flag is True:
                        print('Discarding {0} because it does not contain the '
                              'requested variant {1}.'.format(filepath, ','.join(parameters_dict['variant'])))
                    continue

            if 'run' in parameters_dict and parameters_dict['run']:
                if 'run' in filepath_bids_dictionary:
                    if filepath_bids_dictionary['run'] not in parameters_dict['run']:
                        if debug_flag is True:
                            print('Discarding {0} because its run {1} is not in the list '
                                  'provided {2}'.format(filepath, filepath_bids_dictionary['run'],
                                                        ','.join(parameters_dict['run'])))
                        continue
                else:
                    if debug_flag is True:
                        print('Discarding {0} because it does not contain the '
                              'requested run {1}.'.format(filepath, ','.join(parameters_dict['run'])))
                    continue

            dictionary_key = "_".join(["-".join([k, filepath_bids_dictionary[k]]) for k in
                                       ['pipeline', 'derivative', 'variant', 'acq', 'filetype'] if
                                       k in filepath_bids_dictionary])

            if 'participant' in args.analysis_level:
                dictionary_key += "_" + "-".join(['sub', filepath_bids_dictionary['sub']])

            if 'ses' in filepath_bids_dictionary:
                if parameters_dict['consolidate_sessions'] is True:
                    if parameters_dict['session']:
                        dictionary_key += "_" + "-".join(['ses', '+'.join(parameters_dict['session'])])
                    else:
                        dictionary_key += "_" + "-".join(['ses', 'all'])
                else:
                    dictionary_key += "_" + "-".join(['ses', filepath_bids_dictionary['ses']])

            if 'run' in filepath_bids_dictionary:
                if parameters_dict['consolidate_runs'] is True:
                    if parameters_dict['run']:
                        dictionary_key += "_" + "-".join(['run', '+'.join(parameters_dict['run'])])
                    else:
                        dictionary_key += "_" + "-".join(['run', 'all'])
                else:
                    dictionary_key += "_" + "-".join(['run', filepath_bids_dictionary['run']])

            if parameters_dict['consolidate_tasks'] is True:
                if parameters_dict['task']:
                    dictionary_key += "_" + "-".join(['run', '+'.join(parameters_dict['run'])])
                else:
                    dictionary_key += "_" + "-".join(['run', 'all'])
            else:
                dictionary_key += "_" + "-".join(['task', filepath_bids_dictionary['task']])

            if dictionary_key not in bids_dictionary:
                bids_dictionary[dictionary_key] = []
            bids_dictionary[dictionary_key].append(filepath_bids_dictionary['full_path_name'])
            participant_list.append(filepath_bids_dictionary['participant_label'])

            inclusion_file_count += 1

        if debug_flag is True:
            print("Found {0} files organized into {1} tranches.".format(inclusion_file_count,
                                                                        len(bids_dictionary.keys())))
            print(bids_dictionary)

    # reduce down to just the tranche that the user would like to process now
    if args.participant_ndx:
        participant_index = int(args.participant_ndx)
        if 0 <= int(participant_index) < len(bids_dictionary.keys()):
            participant_key = list(bids_dictionary)[participant_index]
            # make sure to keep it a list
            bids_dictionary = {participant_key: bids_dictionary[participant_key]}
            participant_list = [participant_list[participant_index]]
        else:
            raise ValueError("Participant ndx {0} is out of bounds [0,{1})".format(participant_index,
                                                                                   len(bids_dictionary.keys())))

    if 'test' in args.analysis_level or debug_flag is True:

        # update config file
        if "s3://" not in args.output_dir.lower():
            config_file = os.path.join(args.output_dir, "pyclusterroi_config_{0}.json".format(timestamp))
        else:
            config_file = "pyclusterroi_config_{0}.json".format(timestamp)

        with open(config_file, 'w') as f:
            json.dump(parameters_dict, f)

        # write out the data configuration file
        if args.participant_ndx:
            data_config_file = "pyclusterroi_data_config_pt{0}_{1}.yml".format(args.participant_ndx, timestamp)
        else:
            data_config_file = "pyclusterroi_data_list_{0}.yml".format(timestamp)

        if "s3://" not in args.output_dir.lower():
            data_config_file = os.path.join(args.output_dir, data_config_file)
        else:
            data_config_file = os.path.join(os.getcwd(), data_config_file)

        with open(data_config_file, 'w') as out_stream:
            json.dump(bids_dictionary, out_stream, indent=4)

        if 'test' in args.analysis_level:

            print(
                'System configuration test complete, pipeline configuration written to {0} and data list written '
                'to {1}'.format(config_file, data_config_file))

            return 0

    if os.path.basename(args.output_dir.lstrip('/')) is not 'derivatives':
        output_prefix = os.path.join(args.output_dir, 'derivatives/pyclusterroi')
    else:
        output_prefix = os.path.join(args.output_dir, 'pyclusterroi')

    # if group level analysis, create a directory for the outputs and a .json file mapping the group to participants
    if 'group' in args.analysis_level:

        output_prefix = os.path.join(output_prefix, '-'.join(['grp', parameters_dict['group_name']]))
        os.makedirs(output_prefix, exist_ok=True)

        with open(os.path.join(output_prefix, 'group_definition.json'), 'w') as out_stream:
            json.dump(participant_list, out_stream, indent=4)

    for tranche_key, tranche_paths in bids_dictionary.items():

        tranche_bids = {value.split('-')[0]: value.split('-')[1] for value in tranche_key.split('_')}

        # replace derivative and update variant
        if 'variant' not in tranche_bids:
            tranche_bids['variant'] = ''

        tranche_bids['derivative'] = 'roi'

        if 'group' in args.analysis_level:
            tranche_bids['variant'] += parameters_dict['group_level_method'] + parameters_dict[
                'similarity_metric'] + '{k}'

            tranche_output_prefix = os.path.join(output_prefix,
                                                 'func',
                                                 '-'.join(['grp', parameters_dict['group_name']]) +
                                                 '_'.join(["-".join([bids_key, tranche_bids[bids_key]]) for bids_key in
                                                           ['task', 'run', 'acq'] if bids_key in tranche_bids] +
                                                          [tranche_bids['filetype']] +
                                                          ["-".join([bids_key, tranche_bids[bids_key]]) for bids_key in
                                                           ['variant'] if bids_key in tranche_bids] +
                                                          [tranche_bids['derivative']]))

            print('Performing {0} group level parcellation for:\n{1}\nusing file(s):\n{2}.'.format(
                parameters_dict["group_level_method"], tranche_output_prefix, '\n'.join(tranche_paths)))

        elif 'participant' in args.analysis_level:
            tranche_bids['variant'] += parameters_dict['similarity_metric'] + '{k}'

            tranche_output_prefix = os.path.join(output_prefix,
                                                 "/".join(["-".join([bids_key, tranche_bids[bids_key]]) for bids_key in
                                                           ['sub', 'ses'] if
                                                           bids_key in tranche_bids]),
                                                 'func',
                                                 "_".join(["-".join([bids_key, tranche_bids[bids_key]]) for bids_key in
                                                           ['sub', 'ses', 'task', 'run', 'acq'] if
                                                           bids_key in tranche_bids] + [tranche_bids['filetype']] + [
                                                              "-".join([bids_key, tranche_bids[bids_key]]) for bids_key
                                                              in ['variant'] if bids_key in tranche_bids] + [
                                                              tranche_bids['derivative']]))

            print('Performing participant level parcellation for:\n{0}\nusing file(s):\n{1}.'.format(
                tranche_output_prefix, '\n'.join(tranche_paths)))

        os.makedirs(os.path.dirname(tranche_output_prefix), exist_ok=True)

        out_files = pc.make_parcellation(parameters_dict['number_of_clusters'],
                                         parameters_dict['similarity_metric'],
                                         tranche_paths,
                                         parameters_dict['mask'],
                                         tranche_output_prefix,
                                         thresh=parameters_dict['thresh'],
                                         num_threads=parameters_dict["number_of_threads"])

        if debug_flag is True:
            print("Files generated:\n{0}".format("\n".join(out_files)))

    return 0

    #
    #     with MyPool(num_pt_at_once) as pool:
    #         thread_results = []
    #         for participant in [1, 2, 3]:
    #             participant_input_data = os.path.dirname(
    #                 __file__) + "/test_data/sub-{0}/func/sub-{0}_task-rest_bold.nii.gz".format(participant)
    #
    #             participant_output_prefix = os.path.dirname(
    #                 __file__) + "/test_data/derivatives/pyclusterroi/" + "sub-{0}/func/sub-{0}_task-rest_bold_".format(
    #                 participant) + "variant-tcorr{k}_roi.nii.gz"
    #
    #             os.makedirs(os.path.dirname(participant_output_prefix), exist_ok=True)
    #
    #             thread_results.append(pool.apply_async(pc.make_parcellation, (100, 'tcorr', participant_input_data,
    #                                                                           gm_mask_file, participant_output_prefix,
    #                                                                           0.5, num_threads)))
    #
    #         participant_parcellation_data = []
    #         for thread_result in thread_results:
    #             participant_parcellation_data.append(thread_result.get()[0])
    #
    #     print(participant_parcellation_data)
    #
    #     group_output_prefix = os.path.dirname(
    #         __file__) + "/test_data/derivatives/pyclusterroi/group_task-rest_bold_variant-tcorr{k}_roi.nii.gz"
    #
    #     out_files = pc.make_parcellation(range(10, 100, 10), 'cluster', participant_parcellation_data, gm_mask_file,
    #                                      group_output_prefix, thresh=0.5, num_threads=1)
    #
    #     self.assertTrue(out_files)