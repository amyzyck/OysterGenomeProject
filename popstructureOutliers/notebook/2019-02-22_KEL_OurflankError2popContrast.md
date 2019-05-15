# 20190222 KEL

Bodie kept getting an "pi0" error with a 2 population contrast in OutFlank with the 50K SNP set. We traced the problem to the q-value function. I have seen this problem before with q-value when the p-values either follow the null distribution perfectly or are anit-conservative:

Example of error:
```
> qvalue(p, fdr.level = 0.05, pi0.method="bootstrap")

Error in quantile.default(pi0, prob = 0.1) : 
missing values and NaN's not allowed if 'na.rm' is FALSE
In addition: Warning messages:
1: In min(p) : no non-missing arguments to min; returning Inf
2: In max(p) : no non-missing arguments to max; returning -Inf
3: In min(p) : no non-missing arguments to min; returning Inf
4: In max(p) : no non-missing arguments to max; returning -Inf
```

Note: also get above error when `p=NULL`.

We solved the problem by setting `pi0=1`, which is equivalent to the Benjamini-Hotchberg procedure in the q-value manual and worked on Bodie's computer:

```
qvalue(p, fdr.level = 0.05, pi0=1)
```
