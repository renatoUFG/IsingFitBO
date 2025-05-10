# GNU GENERAL PUBLIC LICENSE - Version 2, June 1991
# Copyright (C) Claudia van Borkulo
# Modified by Renato Rodrigues Silva (2025) to include Bayesian Optimization.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
print.IsingFitBO <-
function(x, ...)
{
  cat("Estimated network:\n")
  
  print(round(x$weiadj,2))
  
  cat("\n\nEstimated Thresholds:\n")
  
  print(x$thresholds)  
}
