\name{madras}
\docType{data}
\alias{madras}
\title{Madras Longitudinal Schizophrenia Study data}
\description{
Data from a study of the natural history of schizophrenia.  First-episode patients had disease symptoms recorded monthly for up to 12 months following initial hospitalization.
}
\usage{data(madras)}
\format{A matrix containing 8 variables (columns) and 922 observations (rows) measured on a total of 86 patients.

\item{id:}{
a vector that identifies the clusters which correspond to the binary response vector given by \code{y}.
}
\item{y:}{
the binary response vector for presence or absence of schizophrenic symptoms
}
\item{month:}{
symptom measurements were collected over the course of 12 months (0,...,11) following the initial hospitalization for schizophrenia
}
\item{age:}{
a binary indicator for age at onset (1 = age < 20 years old)
}
\item{gender:}{
0 = male, 1 = female
}
\item{monthXage:}{
an interaction term (cross product) of month and age
}
\item{monthXgender:}{
an interaction term (cross product) of month and gender
}
\item{initial:}{
an indicator variable for \code{month = 1} allowing \code{alpha1} to be used in lag-2 transition models for both the second-order model and the initial state
}

}
\references{
Thara, R., Henrietta, M., Joseph, A., Rajkumar, S., and Eaton, W. (1994). Ten year course of schizophrenia - The Madras Longitudinal Study. \emph{Acta Psychiatrica Scandinavica} {\bf 90}, 329-336.
}

\keyword{datasets}

