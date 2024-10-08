\name{dftse}
\alias{dftse}
\title{Remove irrelevant frequencies}
\description{
Remove irrelevant frequencies
}
\usage{
dftse(x, low_freq = NULL, high_freq = NULL)
}
\arguments{
\item{x}{
Vector, \code{data.frame}, \code{matrix} or any similar 1D/2D object containing
values for filtering.
}
\item{low_freq}{
Number indicating the lowest period of oscillation as fractions of \eqn{\pi}.
If \code{low_freq > 1}, indicating that the direct frequency of the data is
provided, this is transformed internally into \code{2 / high_freq}. The default
is NULL, meaning that the \code{ifelse(freq > 1, trunc(freq * 1.5), 2)} will be
used.
}
\item{high_freq}{
Number indicating the highest period of oscillation as radians of \eqn{\pi}. If
\code{high_freq > 1}, indicating that the direct frequency of the data is
provided, this is transformed internally into \code{2 / low_freq}. The default
is NULL, meaning that the \code{trunc(freq * 8)} will be used.
}
}
\details{
This is a pure R implementation of removing the irrelevant frequencies. First,
DFT is applied on the data and this result is filtered according to
\code{low_freq} and \code{high_freq}. Finally, an inverse DFT is performed on
these relevant frequencies. Both \code{low_freq} and \code{high_freq} must be
either between 0 and 1, meaning that they are frequencies of the period as
radians, or both >1, indicating that both are starting and ending periods of the
cycle.

\code{low_freq} and \code{high_freq} are used for keeping the relevant
frequencies. These are meant to be the ones inside the range
\eqn{[ low \_ freq, high \_ freq ]}. Therefore, values outside this range are
removed.

For 2-dimensional objects \code{x}, this transformation is applied per column.
}

\value{
Filtered object with length/dimensions same with the input x. Note that for
inputs with dimensions (e.g. \code{matrix}, \code{data.frame}) a \code{matrix}
object will be returned.
}

\references{
Corbae, D., Ouliaris, S., & Phillips, P. (2002), Band Spectral Regression with
Trending-Data. \emph{Econometrica} \bold{70}(3), pp. 1067-1109.

Corbae, D. & Ouliaris, S. (2006), Extracting Cycles from Nonstationary Data,
in Corbae D., Durlauf S.N., & Hansen B.E. (eds.). \emph{Econometric Theory and
Practice: Frontiers of Analysis and Applied Research}. Cambridge: Cambridge
University Press, pp. 167–177. \doi{10.1017/CBO9781139164863.008}.

Shaw, E.S. (1947), Burns and Mitchell on Business Cycles. \emph{Journal of
Political Economy}, \bold{55}(4): pp. 281-298. \doi{10.1086/256533}.
}

\seealso{
\code{\link[corbouli]{corbae_ouliaris}}
}

\examples{
# Apply on ts object
data(USgdp)
res <- dftse(USgdp, low_freq = 0.0625, high_freq = 0.3333)
head(res)

# Apply on vector
res <- dftse(c(USgdp), low_freq = 0.0625, high_freq = 0.3333)
head(res)

# Apply on matrix per column
mat <- matrix(USgdp, ncol = 4)
res <- dftse(mat, low_freq = 0.0625, high_freq = 0.3333)
head(res)

# Apply on data.frame per column
dfmat <- as.data.frame(mat)
res <- dftse(dfmat, low_freq = 0.0625, high_freq = 0.3333)
head(res)
}
