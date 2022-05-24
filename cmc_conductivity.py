import math
import statistics

def determine_cmc(concentrations, conductivities):
    lower_rsqaured_optimized = 0
    lower_linearregion_index_inclusive = 2
    for i in range(3, len(concentrations)):
        lower_rsquared = statistics.correlation(concentrations[0:i], conductivities[0:i])
        if lower_rsquared >= lower_rsqaured_optimized:
            lower_rsqaured_optimized = lower_rsquared
            lower_linearregion_index_inclusive = i - 1
    higher_rsqaured_optimized = 0
    higher_linearregion_index_inclusive = lower_linearregion_index_inclusive + 1
    for i in range(3, len(concentrations)):
        higher_rsquared = statistics.correlation(concentrations[len(concentrations)-1:len(concentrations)-(i+1)])
        if higher_rsquared >=  higher_rsqaured_optimized:
            higher_rsqaured_optimized = higher_rsquared
            higher_linearregion_index_inclusive = len(concentrations) - i
    lower_slope, lower_intercept = statistics.linear_regression(concentrations[0:lower_linearregion_index_inclusive+1], \
         conductivities[0:lower_linearregion_index_inclusive+1])
    higher_slope, higher_intercept = statistics.linear_regression(concentrations[len(concentrations)-1:higher_linearregion_index_inclusive-1], \
         conductivities[len(concentrations)-1:higher_linearregion_index_inclusive-1])
    cmc = (higher_intercept - lower_intercept) / (lower_slope - higher_slope)
    return cmc, lower_slope, lower_intercept, lower_rsqaured_optimized, higher_slope, higher_intercept, higher_rsqaured_optimized
        







