---
layout: page
title: pyClusterROI
tagline: Data driven construction of ROIs from fMRI data. 
---
{% include JB/setup %}

pyClusterROI is a set of python scripts that implement the spatially
constrained clustering methods described in <pubmed>21769991</pubmed>, which is a
method for generating ROIs based on the group-level clustering of functional
MRI data. A spatial constraint is imposed to ensure that the resulting ROIs are
spatially coherent, i.e. the voxels in the resulting ROIs are connected. Using
this package, clustering can be performed based on either the temporal
correlation between voxel time courses, the spatial correlation between whole
brain functional connectivity maps generated from each voxel time course, or a
random approach. Group level clustering can be achieved by either clustering
the average of individual connectivity maps, or a 2-level approach in which
single subject data is clustered, combined, and then submitted to another
clustering. These methods require the specification of the number of clusters
(ROIs) that the user would like to generate. 

Based on the evaluations preformed in our paper, 2-level group clustering using
temporal correlation performed better than other methods. But the user might
choose a different approach based on the specific analysis that is being
performed. The group-mean approach requires much less computation, so it might
be more appropriate for large datasets. The user might prefer spatial
correlation if they specifically want to optimize for the homogeneity of FC
maps generated from withen-ROI clusters. Although, evidence from the paper
suggests that temporal correlation does a better job of optimizing the
homogeneity of FC maps then does spatial correlation. Additionally the number
of clusters generated must be determined by the type of analysis to be
performed. If the desire is to reduce the dimensionality to a low number while
preserving functional homogeneity and interpretability, then a clustering in
the range of 150 to 200 might be optimal. On the other hand, if the desire is
to provide a modest amount of dimensionality reduction, but still preserve
information present at the voxel scale, 600 - 1000 ROIs might be more
appropriate. 

More information about the clustering approach can be found in
<pubmed>21769991</pubmed> as well as a poster (<nitrc group="cluster_roi"
doc="1036" />).


Read [Jekyll Quick Start](http://jekyllbootstrap.com/usage/jekyll-quick-start.html)

Complete usage and documentation available at: [Jekyll Bootstrap](http://jekyllbootstrap.com)

## Update Author Attributes

In `_config.yml` remember to specify your own data:
    
    title : My Blog =)
    
    author :
      name : Name Lastname
      email : blah@email.test
      github : username
      twitter : username

The theme should reference these variables whenever needed.
    
## Sample Posts

This blog contains sample posts which help stage pages and blog data.
When you don't need the samples anymore just delete the `_posts/core-samples` folder.

    $ rm -rf _posts/core-samples

Here's a sample "posts list".

<ul class="posts">
  {% for post in site.posts %}
    <li><span>{{ post.date | date_to_string }}</span> &raquo; <a href="{{ BASE_PATH }}{{ post.url }}">{{ post.title }}</a></li>
  {% endfor %}
</ul>

## To-Do

This theme is still unfinished. If you'd like to be added as a contributor, [please fork](http://github.com/plusjade/jekyll-bootstrap)!
We need to clean up the themes, make theme usage guides with theme-specific markup examples.


