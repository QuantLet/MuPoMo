
[<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/banner.png" width="888" alt="Visit QuantNet">](http://quantlet.de/)

## [<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/qloqo.png" alt="Visit QuantNet">](http://quantlet.de/) **MuPoMo_normalization** [<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/QN2.png" width="60" alt="Visit QuantNet 2.0">](http://quantlet.de/)

```yaml

Name of QuantLet : MuPoMo_normalization

Published in : MuPoMo

Description : 'Normalizes optimal shape variation parameters theta estimated from
MuPoMo_optimization, and will be called in MuPoMo_main_multipop.'

Keywords : time series, demography, mortality, population, normalization

See also : 'MuPoMo_data, MuPoMo_optimization, MuPoMo_referencecurve, MuPoMo_main_twopop,
MuPoMo_main_multipop, MuPoMo_bootstrap'

Author : Lei Fang

Submitted : Mon, June 13 2016 by Lei Fang

```


### R Code:
```r
normalization  = function(theta) {
    theta.temp = colSums(theta)
    for (i in 1:loop.31) {
        nam10  = paste("normal.theta", names.31[i], 1, sep = ".")
        assign(nam10, loop.31 * theta[i, 1]/theta.temp[1])
        nam11  = paste("normal.theta", names.31[i], 3, sep = ".")
        assign(nam11, loop.31 * theta[i, 3]/theta.temp[3])
        nam12  = paste("normal.theta", names.31[i], 2, sep = ".")
        assign(nam12, theta[i, 2] - theta.temp[2]/loop.31)
        nam14  = paste("normal.theta", names.31[i], sep = ".")
        assign(nam14, c(eval(parse(text = paste("normal.theta", names.31[i], 1, sep = "."))), eval(parse(text = paste("normal.theta", 
            names.31[i], 2, sep = "."))), eval(parse(text = paste("normal.theta", names.31[i], 3, sep = ".")))))
    }
    ll = list()
    for (i in 1:loop.31) {
        ll[[i]]  = c(eval(parse(text = paste("normal.theta", names.31[i], sep = "."))))
    }
    normal.theta = do.call(rbind, ll)
    return(normal.theta)
}

```
