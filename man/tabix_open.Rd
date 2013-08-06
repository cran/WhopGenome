\name{tabix_open}
\alias{tabix_open}
\title{
Open Tabix-indexed file for subsequent access with other tabix_ methods
}
\description{
Open Tabix-indexed file for subsequent access with other tabix_ methods
}
\usage{
tabix_open(filename)
}
\arguments{
  \item{filename}{String, name of tabix-indexed file to open}
}
\details{
As filename, specify the data file, not the index file ending in .tbi!
}
\value{
Tabix file handle
}
\author{
Ulrich Wittelsbuerger
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
tabix_close, tabix_read
}
