\name{tabix_getregion}
\alias{tabix_getregion}
\title{
Return the currently selected region of the given tabix file.
}
\description{
Return the currently selected region of the given tabix file. The resulting value does
not reflect the current read position inside that region, i.e. you cannot infer whether
there are any lines left for reading from that region.
}
\usage{
tabix_getregion( tabfh )
}
\arguments{
  \item{tabfh}{Tabix handle, once returned by tabix_open}
}
\details{
Use .Call("tabix_getRegion", tabfh ) to eliminate the slight overhead of using the R wrapper function.
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
