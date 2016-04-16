try_again_message <- function (times, code, message_times = TRUE) 
{
    ## The code of try_again, in testthat v.1.0, with a message for times
    ## run. See my entry in https://github.com/hadley/testthat/issues/451
    init_times <- times
    while (times > 0) {
        e <- tryCatch(withCallingHandlers({
            code
            NULL
        }, warning = function(e) {
            if (identical(e$message, "restarting interrupted promise evaluation")) {
                invokeRestart("muffleWarning")
            }
        }), expectation_failure = function(e) {
            e
        }, error = function(e) {
            e
        })
        if (is.null(e)) {
            if(message_times)
                message("\n times run: ", init_times - times + 1)
            ## is something funny?
            cat(paste("\n cat sentinel time run: ", init_times - times + 1, "\n"))
            return(invisible(TRUE))
        }
        times <- times - 1L
    }
    stop(e)
}
