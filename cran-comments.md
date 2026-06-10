## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.
* To fix the issues that were present in our latest CRAN release we:
    - Rewrote the title of the references in the description in quotes.
    - Added a @return tag to the `print.community_object()` function and updated the RD-tags of print.community_object.Rd: \value.
    - Updated all functions that use the `print()` function to display messages to prefer functions that can be more easily suppressed (e.g. `message()`, `warning()`, or `stop()`).
    - We adjusted the vignette to reset our graphical parameters (based on our use of `par()`).
    - Listed Anton Pervukhin, Florian Rasche, Henner Sudek, Marcel Martin, and Yuanyue Li as copyright holders (cph) and contributors (ctb).
        - We used their code based on there solution on GitHub and this should give them the credit for our use of their solutions.
* I was told to update the print function in the R/community_object.R file, however, the only place we used the print function is the `print.community_object()` function. This function is a S3 generic print function and per the cran-cookbook should be accepted.