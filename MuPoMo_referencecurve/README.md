
[<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/banner.png" width="888" alt="Visit QuantNet">](http://quantlet.de/)

## [<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/qloqo.png" alt="Visit QuantNet">](http://quantlet.de/) **MuPoMo_referencecurve** [<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/QN2.png" width="60" alt="Visit QuantNet 2.0">](http://quantlet.de/)

```yaml

Name of QuantLet : MuPoMo_referencecurve

Published in : MuPoMo

Description : 'Regenerates reference curve (common trend) based on normalized optimal theta
parameters and smoothed original kt, and will be called in MuPoMo_main_multipop.'

Keywords : time series, demography, mortality, population, normalization

See also : 'MuPoMo_data, MuPoMo_optimization, MuPoMo_normalization, MuPoMo_main_twopop,
MuPoMo_main_multipop, MuPoMo_bootstrap'

Author : Lei Fang

Submitted : Mon, June 13 2016 by Lei Fang

```


### R Code:
```r
referencecurve  = function(theta2, theta3, kt, kt.null) {
    t           = time(kt)
    sm.t        = theta3 * t + theta2  # time adjustment
    t.reference = Sweden$year
    ## common grid for kt and shifted kt.reference
    tmin = max(min(t.reference), min(sm.t))
    tmax = min(max(t.reference), max(sm.t))
    i0   = which(t.reference >= tmin & t.reference <= tmax)
    if (length(i0) > 0) {
        t0  = t.reference[i0]
        d   = data.frame(kt, sm.t)
        # sm: sm.regression with optimal smoothing parameter
        h.optimal6 = h.select(sm.t, kt)
        sm  = sm.regression(sm.t, kt, h = h.optimal6, eval.points = t0, model = "none", poly.index = 1, display = "none")
        mu1 = sm$estimate
        mu  = ts(mu1, start = t0[1], frequency = 1)
    } else {
        mu  = kt.null
    }
    return(mu)
}

```
