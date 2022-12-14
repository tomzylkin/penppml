# penppml 0.2.1.900

* When post = FALSE, there was a mistake in the last version where a non-existent
* vector of coefficients was being stored. This issue is now corrected.

# penppml 0.2.0.900

* bootstrap function added
* gamma_val option added for plugin Lasso
* Two possible collinearity checks introduced

# penppml 0.1.1.900

* When standardize option = TRUE and selected penalty is ridge, there was a mistake in last
* version when rescaling the coefficients. The rescaling is now performed correctly. 
* collapse::fhdwithin has replaced lfe::demeanlist in internal hdfeppml_int()-function.
* This has no practical implications for users.

# penppml 0.1.0.900

* Added a `NEWS.md` file to track changes to the package.
