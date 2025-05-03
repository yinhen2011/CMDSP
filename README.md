# CMDSP
[![DOI](https://zenodo.org/badge/951099404.svg)](https://doi.org/10.5281/zenodo.15049533)
The source code is for the paper:

CMDSP: A completed multi-domain shrinkage pattern for texture image classification,

by Bin Li, Xiaochun Xu

PeerJ Computer Science, 2025,

libin126email@126.com；xuxiaochun0303@126.com

version 1.0 (2025.5)

==================================================

Running Environment: Windows 11 pro, Matlab R2024a

==================================================

Reproduce the experimental results for the paper:

You also can download the Outex_TC_00010 dataset from http://www.outex.oulu.fi/index.php?page=classification;
Or you can download the Outex_TC_00010 dataset from Baidu cloud disk link(permanent)-- https://pan.baidu.com/s/1QikrYrFDDOHeoZ31cXscTQ?pwd=tdnk extraction code: tdnk 

In the downloaded file 'Outex_TC_00010', the sub-file 'images' includes all the training and test images, and the sub-file '000' incudes the documents specifying the split of the training and test sets.

UIUC dataset can download from https://slazebni.cs.illinois.edu//.

UMD dataset can download from https://users.umiacs.umd.edu/~fer/High-resolution-data-base/hr_database.htm

Run DomeCodesForOutex10_CMDSP.m to reproduce the reported results. We have added detailed comments to the source code and uploaded the pre-generated multi-scale thresholds for reproducibility. 

Our proposed method is implemented in CMDSP.m, which is a function that implements the LCP and LSP in the paper. 

Cth_3.m and Mth18_3.m are used for generate multi-scale thresholds for LCP and LSP. Before running DomeCodesForOutex10_CMDSP.m, it is need to run Cth_3.m and Mth18_3.m. 

FindUniform.m is used to determine whether the binary pattern is uniform or not, and if it is, it is set to 1, otherwise it is set to 0. It is from the Uniformity Optimization Mechanism (UOM) in the paper.

## Citation

Please consider citing our paper, if you found our work interesting and useful.
```
@article{CMDSP,
  title={Enhanced Texture Classification through a Completed Multi-Domain Shrinkage Pattern},
  author={Bin Li, Xiaochun Xu, and Q.M.Jonathan Wu},
  journal={The Visual computer},
  url={https://github.com/yinhen2011/CMDSP/},
  year={2025}
}
