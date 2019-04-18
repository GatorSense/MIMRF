# MIMRF:
**Multi-Resolution Multi-Modal Sensor Fusion For Remote Sensing Data With Label Uncertainty**

_Xiaoxiao Du and Alina Zare_

If you use this code, cite it: Xiaoxiao Du & Alina Zare. (2019, April 12). GatorSense/MIMRF: Initial Release (Version v1.0). Zenodo. http://doi.org/10.5281/zenodo.2638382 [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2638382.svg)](https://doi.org/10.5281/zenodo.2638382)

[[`arXiv`](https://arxiv.org/abs/1805.00930)] [[`BibTeX`](#CitingMIMRF)]


In this repository, we provide the papers and code for the Multiple Instance Multi-Resolution Fusion (MIMRF) Algorithm.

## Installation Prerequisites

This code uses MATLAB Statistics and Machine Learning Toolbox,
MATLAB Optimization Toolbox and MATLAB Parallel Computing Toolbox.

## Demo

Run `demo_main.m` in MATLAB.

## Main Functions

The MIMRF Algorithm runs using the following function:

```[measure,initialMeasure, Analysis] = learnCIMeasure_minmax_multires(Bags, Labels, Parameters,trueInitMeasure)```


## Inputs

#The *TrainBags* input is a 1xNumTrainBags cell. Inside each cell, NumPntsInBag x nSources cell. Inside each cell, the "collection" of all possible combinations generated from the multi-resolution data set. Details please see Section 3 of the MIMRF paper.

#The *TrainLabels* input is a 1xNumTrainBags double vector that takes values of "1" and "0" for two-class classfication problems -- Training labels for each bag.


## Parameters
The parameters can be set in the following function:

```[Parameters] = learnCIMeasureParams();```

The parameters is a MATLAB structure with the following fields:
1. nPop: size of population
2. sigma: sigma of Gaussians in fitness function
3. maxIterations: maximum number of iterations
4. fitnessThresh: fitness threshold
5. eta: percentage of time to make small-scale mutation
6. sampleVar: variance around sample mean
7. mean: mean of CI in fitness function. This value is always set to 1 (or very close to 1) if the positive label is "1".
8. analysis: if ="1", save all intermediate results

*Parameters can be modified by users in [Parameters] = learnCIMeasureParams() function.*

## Inventory

* Note: some of the util functions were also used in our MICI algorithm. Check out the repo here: [[`MICI Repository`](https://github.com/GatorSense/MICI)]

* Note: the CI-QP (CI fusion using quadratic programming) approach was also implemented and available in our MICI repository.

```
https://github.com/GatorSense/MIMRF

└── root dir
    ├── demo_main.m   //Run this. Main demo file.
    ├── demo_MultiRes_data_MU.mat //Demo multi-resolution dataset
    ├── generateSimData_MU.m //Generates synthetic five-source multi-resolution data
    ├── learnCIMeasureParams.m  //parameters function
    ├── MIMRF_Paper.pdf  //related publication
    ├── learnCIMeasure_minmax_multires.m //MIMRF fusion learning function (learns a fuzzy measure from bag-level data and labels)
    ├── evalFitness_minmax_multires.m //MIMRF fusion fitness evaluation
    ├── computeTestMap.m //MIMRF fusion stage (after learning the optimal fuzzy measure g*, compute CI fusion results)
    └── util  //utility functions
        ├── ChoquetIntegral_g_MultiSources.m  //compute CI for multiple sources
        ├── computeci.m    //compute CI fusion output
        ├── evalInterval.m    //evaluate valid intervals of a fuzzy measure
        ├── ismember_findrow_mex.c  //find row index if vector A is part of a row in vector B.   *Need to run "mex ismember_findrow_mex.c"*
        ├── ismember_findrow_mex_my.m  // find row index if vector A is part of a row in vector B (uses above c code).
        ├── share.h  //global variable header to be used in computeci.c
        ├── invcdf_TruncatedGaussian.m //compute inverse cdf for Truncated Gaussian
        ├── rowcol.m //compute row and column index given image index
        ├── sampleMeasure.m //sample new measures
        ├── sampleMeasure_Above.m  //sampling a new measure from top-down.
        ├── sampleMeasure_Bottom.m  //sampling a new measure from bottom-up.
        └── sampleMultinomial_mat.m  //sample from a multinomial distribution.

```

## License

This source code is licensed under the license found in the [`LICENSE`](LICENSE) file in the root directory of this source tree.

This product is Copyright (c) 2018 X. Du and A. Zare. All rights reserved.

## <a name="CitingMIMRF"></a>Citing MIMRF

If you use the MIMRF multi-resolution fusion algorithm, please cite the following reference using the following BibTeX entries.
```
@article{du2018multi,
  title={Multi-Resolution Multi-Modal Sensor Fusion For Remote Sensing Data With Label Uncertainty},
  author={Du, Xiaoxiao and Zare, Alina},
  journal={arXiv preprint arXiv:1805.00930},
  year={2018}
}
```

## <a name="Related Work"></a>Related Work

Also check out our MICI (Multiple Instance Choquet Integral) algorithm for classifier fusion and regression!

[[`IEEEXplore (MICI Classifier Fusion and Regression paper)`](https://ieeexplore.ieee.org/document/8528500)]

[[`GitHub Code Repository`](https://github.com/GatorSense/MICI)]
