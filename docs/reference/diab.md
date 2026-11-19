# Diabetes Progression

These data consist of observations on 442 patients, with the response of
interest being a quantitative measure of disease progression one year
after baseline.

There are ten baseline variables: age, sex, body-mass index (`bmi`),
average blood pressure (`map`) and six blood serum measurements.

## Usage

``` r
data("diab")
```

## Format

A data frame with 442 observations on the following 11 variables.

- `prog`: disease progression, a numeric vector

- `age`: age, a numeric vector

- `sex`: integer, a numeric vector

- `bmi`: body mass index, a numeric vector

- `map`: mean arterial blood pressure, a numeric vector

- `tc`: blood serum TC, a numeric vector

- `ldl`: blood serum low-density lipoprotein ("bad cholesterol"), a
  numeric vector

- `hdl`: blood serum high-density lipoprotein ("good cholesterol"), a
  numeric vector

- `tch`: blood serum TCH, a numeric vector

- `ltg`: blood serum lamotrigine, a numeric vector

- `glu`: blood serum glucose, a numeric vector

## Source

The dataset was taken from the web site for Efron & Hastie (2021),
*Computer Age Statistical Inference*,
<https://hastie.su.domains/CASI_files/DATA/diabetes.csv>.

## Details

Efron & Hastie describe their analysis using the centered predictor
variables standardized to unit L2 norm.
[`ridge()`](https://friendly.github.io/genridge/reference/ridge.md) does
not (yet) provide this scaling.

## References

Efron, B., Hastie, T., Johnstone, I., & Tibshirani, R. (2004). Least
Angle Regression. *The Annals of Statistics*, **32**(2), 407-499.
[doi:10.1214/009053604000000067](https://doi.org/10.1214/009053604000000067)

Efron, B., & Hastie, T. (2021). *Computer Age Statistical Inference,
Student Edition: Algorithms, Evidence, and Data Science*, Cambridge
University Press.
[doi:10.1017/9781108914062](https://doi.org/10.1017/9781108914062)

## Examples

``` r
data(diab)
## maybe str(diab) ; plot(diab) ...
```
