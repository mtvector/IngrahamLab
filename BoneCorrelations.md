BoneCorrelations
================
Matthew Schmitz
March 15, 2018

So to get at the question of what is changing, perhaps we can use the variability of the FACS efficiency in the Ambrosi data, and the inter-individual animal variation to our advantage. I'll try to do this by looking for genes which covary across the samples, expecting to find sets related to cell type which are higher in their different sorted populations.

Load all the processed data from BoneNotebook.Rmd's cache

``` r
load("~/code/IngrahamLab/BoneNotebook_cache/markdown_github/everything.RData")
```

    ## 

Apply WGCNA:

``` r
library(WGCNA)
library(matrixStats)
mat <- log(ambrosiMatNorm+1,2)
sum(rowSds(mat)>1)
```

    ## [1] 12939

``` r
softThresh <- WGCNA::pickSoftThreshold(t(mat[rowSds(mat)>1,]))
```

    ## Warning: executing %dopar% sequentially: no parallel backend registered

    ##    Power SFT.R.sq  slope truncated.R.sq mean.k. median.k. max.k.
    ## 1      1 0.065800  2.220          0.691  3520.0   3530.00   4220
    ## 2      2 0.000717  0.128          0.802  1420.0   1430.00   1980
    ## 3      3 0.033000 -0.705          0.986   707.0    704.00   1120
    ## 4      4 0.057200 -0.851          0.972   402.0    391.00    715
    ## 5      5 0.083600 -0.810          0.916   251.0    239.00    490
    ## 6      6 0.130000 -0.649          0.787   169.0    156.00    355
    ## 7      7 0.451000 -1.210          0.845   120.0    108.00    306
    ## 8      8 0.706000 -1.530          0.897    89.7     76.10    275
    ## 9      9 0.838000 -1.690          0.914    69.4     55.30    252
    ## 10    10 0.897000 -1.780          0.920    55.5     41.30    235
    ## 11    12 0.800000 -1.920          0.764    38.2     24.60    212
    ## 12    14 0.784000 -1.850          0.730    28.4     15.70    196
    ## 13    16 0.766000 -1.730          0.704    22.4     10.60    185
    ## 14    18 0.704000 -1.720          0.631    18.4      7.34    177
    ## 15    20 0.724000 -1.610          0.649    15.7      5.30    170

``` r
#hardThresh <- pickHardThreshold(t(log(ambrosiMatNorm+1,2)))
AmbrosiBlocks <- blockwiseModules(t(mat[rowSds(mat)>1,]),maxBlockSize = 7000,power = softThresh$powerEstimate,deepSplit = 1)
AmbrosiMatLog <- mat
```

``` r
mat <- log(boneMatNorm+1,2)
sum(rowSds(mat)>.6)
```

    ## [1] 11785

``` r
softThresh <- WGCNA::pickSoftThreshold(t(mat[rowSds(mat)>.6,]))
```

    ##    Power SFT.R.sq slope truncated.R.sq mean.k. median.k. max.k.
    ## 1      1   0.0838  2.35          0.630  3940.0    3860.0   5120
    ## 2      2   0.2660 -1.47          0.610  1900.0    1810.0   3020
    ## 3      3   0.7120 -2.01          0.825  1090.0     998.0   2060
    ## 4      4   0.8110 -2.07          0.872   700.0     612.0   1530
    ## 5      5   0.8590 -2.06          0.899   481.0     407.0   1190
    ## 6      6   0.8880 -2.07          0.915   349.0     284.0    965
    ## 7      7   0.9080 -2.05          0.924   264.0     212.0    803
    ## 8      8   0.9260 -2.03          0.936   205.0     162.0    682
    ## 9      9   0.9340 -2.02          0.944   164.0     125.0    589
    ## 10    10   0.9380 -2.00          0.947   134.0      99.5    515
    ## 11    12   0.9490 -1.94          0.952    94.1      66.1    409
    ## 12    14   0.9490 -1.89          0.955    69.6      46.3    335
    ## 13    16   0.9380 -1.84          0.940    53.7      33.8    282
    ## 14    18   0.9710 -1.74          0.971    42.8      25.4    243
    ## 15    20   0.9580 -1.68          0.950    35.1      19.5    212

``` r
softThresh
```

    ## $powerEstimate
    ## [1] 5
    ## 
    ## $fitIndices
    ##    Power   SFT.R.sq     slope truncated.R.sq    mean.k.  median.k.
    ## 1      1 0.08378535  2.348680      0.6304164 3943.60024 3860.43710
    ## 2      2 0.26564020 -1.466177      0.6102762 1903.65466 1805.60894
    ## 3      3 0.71173125 -2.012691      0.8253354 1094.50479  997.70773
    ## 4      4 0.81089126 -2.072865      0.8716966  700.06999  611.93389
    ## 5      5 0.85940142 -2.061669      0.8988047  481.49366  406.98336
    ## 6      6 0.88776362 -2.068984      0.9149800  349.12268  283.67793
    ## 7      7 0.90775013 -2.052106      0.9237866  263.53461  212.12711
    ## 8      8 0.92552484 -2.032886      0.9356032  205.33771  161.77091
    ## 9      9 0.93431731 -2.020511      0.9436547  164.15311  125.42544
    ## 10    10 0.93799454 -2.002806      0.9468082  134.04599   99.50354
    ## 11    12 0.94864813 -1.943666      0.9516518   94.06898   66.13212
    ## 12    14 0.94923959 -1.891236      0.9549300   69.63570   46.29969
    ## 13    16 0.93829331 -1.843852      0.9401849   53.72075   33.78666
    ## 14    18 0.97113294 -1.743935      0.9707190   42.83362   25.37598
    ## 15    20 0.95834011 -1.675896      0.9496271   35.09044   19.54236
    ##       max.k.
    ## 1  5119.7829
    ## 2  3023.6760
    ## 3  2060.4337
    ## 4  1526.1545
    ## 5  1191.5776
    ## 6   965.1514
    ## 7   803.1322
    ## 8   682.2482
    ## 9   589.0830
    ## 10  515.4012
    ## 11  408.5679
    ## 12  335.3436
    ## 13  282.3150
    ## 14  242.9579
    ## 15  212.3544

``` r
#hardThresh <- pickHardThreshold(t(log(boneMatNorm+1,2)))
IngrahamBlocks <- blockwiseModules(t(mat[rowSds(mat)>.6,]),maxBlockSize = 7000,power = softThresh$powerEstimate,deepSplit = 1)
IngrahamMatLog <- mat
```

``` r
save.image(file = "~/code/IngrahamLab/BoneCorrelations_cache/wgcnas.RData")
#load(file = "~/code/IngrahamLab/BoneCorrelations_cache/wgcnas.RData")

par(mfrow=c(2,2))
lapply(1:length(IngrahamBlocks$dendrograms),function(i) plotDendroAndColors(IngrahamBlocks$dendrograms[[i]],colors=IngrahamBlocks$colors[IngrahamBlocks$blockGenes[[i]]]))
```

![](BoneCorrelations_files/figure-markdown_github/Trees-1.png)![](BoneCorrelations_files/figure-markdown_github/Trees-2.png)

    ## [[1]]
    ## [[1]]$mar
    ## [1] 1 5 0 1
    ## 
    ## 
    ## [[2]]
    ## [[2]]$mar
    ## [1] 1 5 0 1

``` r
par(mfrow=c(2,2))
lapply(1:length(AmbrosiBlocks$dendrograms),function(i) plotDendroAndColors(AmbrosiBlocks$dendrograms[[i]],colors=AmbrosiBlocks$colors[AmbrosiBlocks$blockGenes[[i]]]))
```

![](BoneCorrelations_files/figure-markdown_github/Trees-3.png)![](BoneCorrelations_files/figure-markdown_github/Trees-4.png)![](BoneCorrelations_files/figure-markdown_github/Trees-5.png)

    ## [[1]]
    ## [[1]]$mar
    ## [1] 1 5 0 1
    ## 
    ## 
    ## [[2]]
    ## [[2]]$mar
    ## [1] 1 5 0 1
    ## 
    ## 
    ## [[3]]
    ## [[3]]$mar
    ## [1] 1 5 0 1

``` r
print(str(IngrahamBlocks))
```

    ## List of 10
    ##  $ colors        : chr [1:11785] "cyan" "cyan" "turquoise" "turquoise" ...
    ##  $ unmergedColors: chr [1:11785] "darkseagreen4" "darkseagreen4" "turquoise" "pink" ...
    ##  $ MEs           :'data.frame':  8 obs. of  99 variables:
    ##   ..$ MEcoral          : num [1:8] -0.4037 -0.3343 0.6003 -0.1355 0.0415 ...
    ##   ..$ MEorange         : num [1:8] -0.181 -0.262 -0.321 0.411 0.453 ...
    ##   ..$ MEplum3          : num [1:8] -0.478 -0.222 0.135 0.151 0.475 ...
    ##   ..$ MEdarkseagreen4  : num [1:8] -0.434 -0.573 -0.33 0.159 0.306 ...
    ##   ..$ MEhoneydew       : num [1:8] -0.186 -0.344 -0.23 -0.251 0.474 ...
    ##   ..$ MEdarkslateblue  : num [1:8] 0.58966 0.59699 -0.19916 0.00477 -0.29074 ...
    ##   ..$ MEorangered4     : num [1:8] 0.6798 -0.1813 -0.0843 0.5251 -0.2583 ...
    ##   ..$ MEred            : num [1:8] 0.931 -0.152 -0.051 -0.158 -0.129 ...
    ##   ..$ MElightpink4     : num [1:8] 0.511 -0.304 -0.217 0.369 -0.346 ...
    ##   ..$ MEbrown4         : num [1:8] 0.388 -0.5 -0.311 0.358 -0.152 ...
    ##   ..$ MEindianred4     : num [1:8] 0.222 -0.51 0.199 0.324 -0.414 ...
    ##   ..$ MEmaroon         : num [1:8] 0.4787 -0.0322 -0.3394 -0.1419 -0.2675 ...
    ##   ..$ MEdarkgrey       : num [1:8] 0.1129 -0.3743 -0.3755 -0.1107 0.0184 ...
    ##   ..$ MEpink4          : num [1:8] 0.36 -0.366 -0.119 -0.443 -0.132 ...
    ##   ..$ MEsienna4        : num [1:8] 0.517 -0.318 0.346 -0.304 -0.362 ...
    ##   ..$ MEyellow3        : num [1:8] 0.432 0.174 0.285 -0.159 -0.304 ...
    ##   ..$ MEgrey60         : num [1:8] -0.122 0.559 -0.118 -0.279 -0.172 ...
    ##   ..$ MEthistle1       : num [1:8] -0.24 0.39 0.249 0.186 -0.313 ...
    ##   ..$ MElightsteelblue : num [1:8] 0.268 0.198 -0.364 0.241 0.199 ...
    ##   ..$ MEcoral3         : num [1:8] 0.331 -0.366 0.257 0.191 0.245 ...
    ##   ..$ MEskyblue2       : num [1:8] 0.181 0.207 0.156 0.137 0.165 ...
    ##   ..$ MEcoral1         : num [1:8] 0.275 -0.241 -0.395 0.256 -0.332 ...
    ##   ..$ MElightcoral     : num [1:8] -0.1922 -0.2537 -0.3087 -0.0562 -0.1255 ...
    ##   ..$ MEfirebrick4     : num [1:8] -0.281 0.338 -0.253 0.376 -0.521 ...
    ##   ..$ MEorangered1     : num [1:8] -0.332 -0.252 -0.294 0.546 -0.248 ...
    ##   ..$ MEorangered3     : num [1:8] -0.399 -0.136 0.156 0.289 -0.37 ...
    ##   ..$ MEplum           : num [1:8] 0.516 -0.22 -0.26 -0.436 -0.15 ...
    ##   ..$ MEskyblue3       : num [1:8] 0.288 0.35 -0.445 -0.374 -0.153 ...
    ##   ..$ MEthistle        : num [1:8] 0.283 -0.219 -0.65 -0.3 -0.124 ...
    ##   ..$ MEviolet         : num [1:8] 0.1888 0.0752 -0.6116 0.1804 -0.5719 ...
    ##   ..$ MEtan            : num [1:8] -0.2009 0.5829 -0.1782 -0.0137 -0.2383 ...
    ##   ..$ MEcoral2         : num [1:8] -0.252 0.274 0.524 -0.301 -0.222 ...
    ##   ..$ MEmagenta4       : num [1:8] -0.328 0.452 0.219 -0.335 -0.307 ...
    ##   ..$ MElavenderblush2 : num [1:8] -0.424 0.242 0.335 0.279 -0.238 ...
    ##   ..$ MEpalevioletred2 : num [1:8] -0.3103 0.5256 0.2112 0.0766 -0.342 ...
    ##   ..$ MElightgreen     : num [1:8] 0.236 0.401 0.323 0.236 -0.397 ...
    ##   ..$ MEwhite          : num [1:8] 0.079 0.148 0.256 0.225 -0.661 ...
    ##   ..$ MEdarkmagenta    : num [1:8] 0.572 -0.158 -0.197 -0.212 -0.238 ...
    ##   ..$ MEblue2          : num [1:8] 0.388 -0.256 0.497 -0.253 -0.32 ...
    ##   ..$ MEdarkviolet     : num [1:8] 0.379 -0.345 0.254 0.325 -0.321 ...
    ##   ..$ MEmediumpurple2  : num [1:8] -0.113 -0.142 0.358 -0.278 0.167 ...
    ##   ..$ MEthistle2       : num [1:8] -0.23 -0.305 0.548 -0.167 -0.223 ...
    ##   ..$ MEbrown          : num [1:8] -0.119 -0.146 -0.178 -0.167 -0.131 ...
    ##   ..$ MEcyan           : num [1:8] -0.155 -0.188 -0.106 0.474 -0.315 ...
    ##   ..$ MEdarkseagreen3  : num [1:8] -0.545 0.225 0.434 0.196 0.422 ...
    ##   ..$ MEivory          : num [1:8] -0.33 0.207 0.345 -0.304 0.549 ...
    ##   ..$ MEpurple         : num [1:8] 0.229 -0.125 0.317 0.158 0.259 ...
    ##   ..$ MEantiquewhite4  : num [1:8] -0.162 -0.142 0.205 0.378 0.216 ...
    ##   ..$ MEdarkorange2    : num [1:8] -0.422 0.233 0.141 0.379 0.249 ...
    ##   ..$ MEhoneydew1      : num [1:8] -0.343 -0.146 -0.197 0.393 0.447 ...
    ##   ..$ MEantiquewhite2  : num [1:8] 0.38 -0.249 -0.183 -0.329 0.366 ...
    ##   ..$ MEroyalblue      : num [1:8] -0.142 -0.21 -0.309 -0.145 0.595 ...
    ##   ..$ MEyellowgreen    : num [1:8] 0.308 0.255 -0.38 0.324 -0.373 ...
    ##   ..$ MEgreenyellow    : num [1:8] 0.102 0.173 -0.558 0.246 0.235 ...
    ##   ..$ MEplum2          : num [1:8] 0.2716 0.0223 -0.4295 0.325 0.235 ...
    ##   ..$ MElightyellow    : num [1:8] -0.236 -0.156 -0.188 -0.206 0.561 ...
    ##   ..$ MEmidnightblue   : num [1:8] 0.195 0.159 0.403 -0.369 0.256 ...
    ##   ..$ MElightsteelblue1: num [1:8] -0.0543 -0.249 0.3389 0.1166 0.3318 ...
    ##   ..$ MEsalmon2        : num [1:8] -0.184 -0.318 0.424 -0.33 0.576 ...
    ##   ..$ MEdarkolivegreen4: num [1:8] -0.267 -0.261 0.152 0.555 0.446 ...
    ##   ..$ MEpalevioletred3 : num [1:8] -0.165 -0.139 -0.213 0.48 0.721 ...
    ##   ..$ MEmediumpurple3  : num [1:8] 0.448 -0.184 -0.241 -0.11 0.741 ...
    ##   ..$ MEyellow         : num [1:8] -0.1435 -0.1577 -0.0607 -0.1082 0.9309 ...
    ##   ..$ MEnavajowhite2   : num [1:8] 0.399 -0.262 -0.183 0.441 0.496 ...
    ##   ..$ MEsalmon         : num [1:8] 0.4761 -0.2555 -0.1531 -0.0601 0.4954 ...
    ##   ..$ MEsaddlebrown    : num [1:8] 0.5 -0.133 -0.208 -0.119 -0.182 ...
    ##   ..$ MEbisque4        : num [1:8] -0.2152 0.0337 0.5242 -0.1949 -0.2997 ...
    ##   ..$ MEmagenta        : num [1:8] -0.1427 -0.0928 -0.0669 -0.0869 -0.1342 ...
    ##   ..$ MEplum1          : num [1:8] -0.1935 0.4309 -0.1083 -0.0851 -0.127 ...
    ##   ..$ MEdarkgreen      : num [1:8] -0.2738 0.0485 0.2278 0.5143 -0.2043 ...
    ##   ..$ MEdarkorange     : num [1:8] -0.154 0.543 -0.219 0.353 -0.02 ...
    ##   ..$ MEnavajowhite1   : num [1:8] 0.2358 0.097 0.0537 0.3816 0.1355 ...
    ##   ..$ MEskyblue        : num [1:8] 0.4025 0.0982 -0.2681 0.5557 -0.2102 ...
    ##   ..$ MEsalmon4        : num [1:8] 0.2591 -0.0745 0.0211 0.2247 -0.2358 ...
    ##   ..$ MEdarkolivegreen : num [1:8] 0.321 0.22 0.163 -0.146 0.142 ...
    ##   ..$ MEturquoise      : num [1:8] 0.1129 0.2046 0.1597 0.1154 0.0702 ...
    ##   ..$ MEmediumorchid   : num [1:8] 0.366 0.524 0.163 -0.209 -0.354 ...
    ##   ..$ MEskyblue4       : num [1:8] 0.4577 0.00388 0.47512 -0.11435 -0.34258 ...
    ##   ..$ MElightpink3     : num [1:8] 0.3822 0.4144 0.3437 -0.3836 0.0327 ...
    ##   ..$ MEdarkturquoise  : num [1:8] 0.154 0.204 0.125 -0.572 0.259 ...
    ##   ..$ MEyellow4        : num [1:8] 0.0615 0.4828 -0.111 -0.3751 0.3327 ...
    ##   ..$ MEskyblue1       : num [1:8] 0.494 0.227 -0.194 -0.564 0.308 ...
    ##   ..$ MEfloralwhite    : num [1:8] 0.428 0.374 -0.201 -0.175 0.503 ...
    ##   ..$ MEthistle3       : num [1:8] 0.425 0.147 0.249 -0.366 0.52 ...
    ##   ..$ MElightcyan1     : num [1:8] -0.147 0.581 -0.145 -0.253 0.636 ...
    ##   ..$ MEmediumpurple1  : num [1:8] -0.254 0.449 -0.196 0.292 0.588 ...
    ##   ..$ MEmediumpurple4  : num [1:8] -0.165 0.718 0.132 0.105 0.188 ...
    ##   ..$ MEpink           : num [1:8] -0.11 0.934 -0.127 -0.166 -0.125 ...
    ##   ..$ MEsienna3        : num [1:8] -0.0805 0.5161 0.4834 0.3037 -0.305 ...
    ##   ..$ MEsteelblue      : num [1:8] -0.116 0.596 0.617 -0.193 -0.194 ...
    ##   ..$ MEblue           : num [1:8] 0.0933 0.4136 0.3546 0.0947 0.1949 ...
    ##   ..$ MEdarkred        : num [1:8] 0.349 0.284 0.189 0.33 0.198 ...
    ##   ..$ MEgreen          : num [1:8] -0.167 -0.114 0.931 -0.1 -0.117 ...
    ##   ..$ MEpaleturquoise  : num [1:8] 0.5345 -0.1234 0.6547 -0.2273 -0.0341 ...
    ##   ..$ MEbrown2         : num [1:8] -0.227 -0.197 0.586 0.627 -0.272 ...
    ##   ..$ MElavenderblush3 : num [1:8] 0.336 -0.249 0.442 0.537 -0.275 ...
    ##   ..$ MEblack          : num [1:8] -0.1235 -0.0784 -0.1061 0.932 -0.1343 ...
    ##   ..$ MElightcyan      : num [1:8] -0.2201 0.5072 0.0186 0.6566 -0.1984 ...
    ##   ..$ MEgrey           : num [1:8] -0.6011 0.0432 0.2377 -0.1129 -0.296 ...
    ##  $ goodSamples   : logi [1:8] TRUE TRUE TRUE TRUE TRUE TRUE ...
    ##  $ goodGenes     : logi [1:11785] TRUE TRUE TRUE TRUE TRUE TRUE ...
    ##  $ dendrograms   :List of 2
    ##   ..$ :List of 7
    ##   .. ..$ merge      : int [1:6751, 1:2] -81 -1497 -1700 -1727 -1874 -4446 -4614 -4746 -4811 -137 ...
    ##   .. ..$ height     : num [1:6751] 0.473 0.475 0.475 0.475 0.475 ...
    ##   .. ..$ order      : int [1:6752] 611 1397 3307 5920 4219 5276 3155 6183 5087 5496 ...
    ##   .. ..$ labels     : NULL
    ##   .. ..$ method     : chr "average"
    ##   .. ..$ call       : language fastcluster::hclust(d = as.dist(dissTom), method = "average")
    ##   .. ..$ dist.method: NULL
    ##   .. ..- attr(*, "class")= chr "hclust"
    ##   ..$ :List of 7
    ##   .. ..$ merge      : int [1:5032, 1:2] -1489 -1599 -1600 -1601 -1612 -1739 -1074 -1087 -1097 -1229 ...
    ##   .. ..$ height     : num [1:5032] 0.423 0.423 0.423 0.423 0.423 ...
    ##   .. ..$ order      : int [1:5033] 3063 4575 2732 3546 1005 2516 683 4181 2308 104 ...
    ##   .. ..$ labels     : NULL
    ##   .. ..$ method     : chr "average"
    ##   .. ..$ call       : language fastcluster::hclust(d = as.dist(dissTom), method = "average")
    ##   .. ..$ dist.method: NULL
    ##   .. ..- attr(*, "class")= chr "hclust"
    ##  $ TOMFiles      : NULL
    ##  $ blockGenes    :List of 2
    ##   ..$ : int [1:6752] 3 4 5 8 9 11 13 16 18 24 ...
    ##   ..$ : int [1:5033] 1 2 6 7 10 12 14 15 17 19 ...
    ##  $ blocks        : num [1:11785] 2 2 1 1 1 2 2 1 1 2 ...
    ##  $ MEsOK         : logi TRUE
    ## NULL

``` r
table(IngrahamBlocks$colors)
```

    ## 
    ##   antiquewhite2   antiquewhite4         bisque4           black 
    ##              38              70              82             254 
    ##            blue           blue2           brown          brown2 
    ##            1035              51             528              52 
    ##          brown4           coral          coral1          coral2 
    ##              83              38              72              69 
    ##          coral3            cyan       darkgreen        darkgrey 
    ##              38             171             127             123 
    ##     darkmagenta  darkolivegreen darkolivegreen4      darkorange 
    ##             103             105              52             121 
    ##     darkorange2         darkred   darkseagreen3   darkseagreen4 
    ##              84             140              41              72 
    ##   darkslateblue   darkturquoise      darkviolet      firebrick4 
    ##              82             125              50              53 
    ##     floralwhite           green     greenyellow            grey 
    ##              85             314             202               1 
    ##          grey60        honeydew       honeydew1      indianred4 
    ##             153              41              73              58 
    ##           ivory  lavenderblush2  lavenderblush3      lightcoral 
    ##              86              42              74              60 
    ##       lightcyan      lightcyan1      lightgreen      lightpink3 
    ##             158              90             148              44 
    ##      lightpink4  lightsteelblue lightsteelblue1     lightyellow 
    ##              74              62              93             147 
    ##         magenta        magenta4          maroon    mediumorchid 
    ##             218              45              76              69 
    ##   mediumpurple1   mediumpurple2   mediumpurple3   mediumpurple4 
    ##              27              64              93              38 
    ##    midnightblue    navajowhite1    navajowhite2          orange 
    ##             163              45              76             122 
    ##      orangered1      orangered3      orangered4   paleturquoise 
    ##              29              66              93             108 
    ##  palevioletred2  palevioletred3            pink           pink4 
    ##              47              78             235              29 
    ##            plum           plum1           plum2           plum3 
    ##              67              94              81              49 
    ##          purple             red       royalblue     saddlebrown 
    ##             203             284             144             115 
    ##          salmon         salmon2         salmon4         sienna3 
    ##             175              47              78             102 
    ##         sienna4         skyblue        skyblue1        skyblue2 
    ##              31             118              67              68 
    ##        skyblue3        skyblue4       steelblue             tan 
    ##              95              37             113             177 
    ##         thistle        thistle1        thistle2        thistle3 
    ##              47              78              79              49 
    ##       turquoise          violet           white          yellow 
    ##            1059             105             120             319 
    ##         yellow3         yellow4     yellowgreen 
    ##              31              68             100
