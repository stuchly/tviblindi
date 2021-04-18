## Printing functions

.msg          <- function(message, endline = TRUE) cat(crayon::bold(crayon::blue(paste0(message, if (endline) '\n' else ''))))
.msg_alt      <- function(message, endline = TRUE) cat(crayon::italic(crayon::bold(paste0(message, if (endline) '\n' else ''))))
.msg_name     <- function(message, endline = TRUE) cat(crayon::bold(crayon::bgWhite(paste0(message, if (endline) '\n' else ''))))
.msg_alt_good <- function(message, endline = TRUE) cat(crayon::green(crayon::italic(crayon::bold(paste0(message, if (endline) '\n' else '')))))
.msg_alt_bad  <- function(message, endline = TRUE) cat(crayon::red(crayon::italic(crayon::bold(paste0(message, if (endline) '\n' else '')))))
.msg_name     <- function(message, endline = TRUE) cat(crayon::bgWhite(crayon::bold(crayon::black(paste0(message, if (endline) '\n' else '')))))
.msg_lite     <- function(message, endline = TRUE) cat(crayon::italic(crayon::black(paste0(message, if (endline) '\n' else ''))))