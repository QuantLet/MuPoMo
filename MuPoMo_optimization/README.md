
[<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/banner.png" width="888" alt="Visit QuantNet">](http://quantlet.de/)

## [<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/qloqo.png" alt="Visit QuantNet">](http://quantlet.de/) **MuPoMo_optimization** [<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/QN2.png" width="60" alt="Visit QuantNet 2.0">](http://quantlet.de/)

```yaml

Name of QuantLet : MuPoMo_optimization

Published in : MuPoMo

Description : 'Optimizes shape variation parameters theta based on kt and reference curve to update
new kt, and will be called in MuPoMo_main_twopop and MuPoMo_main_multipop.'

Keywords : 'time series, demography, mortality, population, nonlinear, optimization, non parametric
smoothing'

See also : 'MuPoMo_data, MuPoMo_normalization, MuPoMo_referencecurve, MuPoMo_main_twopop,
MuPoMo_main_multipop, MuPoMo_bootstrap'

Author : Lei Fang, Juhyun Park

Submitted : Mon, June 13 2016 by Lei Fang

```


### R Code:
```r
optimization = function(theta, kt, kt.reference) {
           t = time(kt)
 t.reference = time(kt.reference)
    
    ### loss function
    loss = function(theta, t, kt, t.reference, kt.reference) {
        ## assume theta[1]>0, theta[3]>0, if not, take absolute values
        theta1 = abs(theta[1])
        theta2 = theta[2]
        theta3 = abs(theta[3])
        
        if (theta[1] < 0 | theta[3] < 0) 
            warning("theta1 or theta3 <0, abs value is used")
        sm.t = (t.reference - theta2)/theta3  # time adjustment
        ## common domain for kt and time-adjusted kt.reference
        tmin = max(min(t), min(sm.t))
        tmax = min(max(t), max(sm.t))
        i0   = which(t >= tmin & t <= tmax)  # index for common domain
        if (length(i0) > 0) {
            t0 = t[i0]
            ## smooth interpolation of shifted kt.reference on common grid
            dref = data.frame(sm.t = sm.t, kt.reference = kt.reference)
            # sm: sm.regression with optimal smoothing parameter
            h.optimal3 = h.select(sm.t, kt.reference)
            sm = sm.regression(sm.t, kt.reference, h = h.optimal3, eval.points = t0, model = "none", poly.index = 1, display = "none")
            mu = theta1 * sm$estimate
            ## mean squared error at common grid points
            mse = mean((kt[i0] - mu)^2)  # mse of kt and the modelled one
        } else {
            mse = 1e+09
        }
        return(mse)
    }
    
    ### parameter estimation with nonlinear optimization
    conv   = 1
    theta0 = theta
    while (conv != 0) {
        # check convergence
        out    = optim(theta0, loss, gr = NULL, t, kt, t.reference, kt.reference, control = list(maxit = 1000))  # optimization
        conv   = out$convergence
        theta0 = out$par
    }
    ### constraint on theta2 (time shift parameter)
    if (out$par[2] >= -200 & out$par[2] <= 200) 
        temp.par = out$par else temp.par = theta
    result = c(temp.par, out$value, out$convergence)
    
    ### check results on graph
    theta  = result[1:3]
    theta1 = abs(theta[1])
    theta2 = theta[2]
    theta3 = abs(theta[3])
    
    sm.t   = (t.reference - theta2)/theta3  # time adjustment
    ## common grid for kt and shifted kt.reference
    tmin   = max(min(t), min(sm.t))
    tmax   = min(max(t), max(sm.t))
    i0     = which(t >= tmin & t <= tmax)
    t0     = t[i0]
    ## shifted curves (kt.hat)
    dref   = data.frame(sm.t = sm.t, kt.reference = kt.reference)
    
    # sm: sm.regression with optimal smoothing parameter
    h.optimal4 = h.select(sm.t, kt.reference)
    sm         = sm.regression(sm.t, kt.reference, h = h.optimal4, eval.points = sm.t, model = "none", poly.index = 1, display = "none")
    mu         = theta1 * sm$estimate
    mu         = ts(mu, start = sm.t[1], frequency = theta3)
    
    # sm: sm.regression with optimal smoothing parameter
    h.optimal5 = h.select(sm.t, kt.reference)
    sm0        = sm.regression(sm.t, kt.reference, h = h.optimal5, eval.points = t0, model = "none", poly.index = 1, display = "none")
    mu0        = theta1 * sm0$estimate
    mu0        = ts(mu0, start = t0[1], frequency = 1)
    
    return(list(result, mu, mu0, tmin, tmax))
}

```
