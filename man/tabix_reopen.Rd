\name{tabix_reopen}
\alias{tabix_reopen}
\title{
Reopen a Tabix-indexed file if the filehandle became invalid.
}
\description{
Reopen a Tabix-indexed file if the filehandle became invalid.
}
\usage{
tabix_reopen( tabfh )
}
\arguments{
  \item{tabfh}{Tabix handle, once returned by tabix_open}
}
\details{
Use .Call("tabix_reopen", tabfh ) to eliminate the slight overhead of using the R wrapper function.
}
\value{
Tabix file handle
}
\author{
Ulrich Wittelsbuerger
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
tabix_open
}
