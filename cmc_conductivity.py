import math
import statistics
import matplotlib.pyplot as plt

TEST_CONCENTRATIONS = [6.003, 6.671, 7.412, \
     8.236, 9.151, 10.167, 11.297, 12.552, 13.947, \
          15.497, 17.219, 19.132, 21.258, 23.620, \
               26.244, 29.160, 32.400, 36.000, 40.000]

TEST_CONDUCTIVITIES = [353, 389, 430, 471, 522, 574, 635, 694, \
    762, 828, 890, 952, 1014, 1079, 1145, 1222, 1304, 1393, 1492]


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
    return cmc, lower_slope, lower_intercept, lower_rsqaured_optimized, higher_slope, higher_intercept, higher_rsqaured_optimized, lower_linearregion_index_inclusive, higher_linearregion_index_inclusive

def plot_cmc(concentrations, conductivities, cmc, lower_slope, lower_intercept, lower_rsqaured_optimized, higher_slope, \
     higher_intercept, higher_rsqaured_optimized, lower_linearregion_index_inclusive, higher_linearregion_index_inclusive):
     plt.scatter(concentrations,conductivities,color="k")
     plt.scatter(concentrations[0:lower_linearregion_index_inclusive+1], conductivities[0:lower_linearregion_index_inclusive+1], "r")
     plt.scatter(concentrations[higher_linearregion_index_inclusive:len(concentrations)], conductivities[higher_linearregion_index_inclusive:len(concentrations)], "g")
     concentration_steprange = concentrations[len(concentrations)-1] - concentrations[0]
     concentration_stepsize = concentration_steprange / 99
     concentration_line_xaxis = []
     concentration_lowerline_yaxis = []
     concentration_higherline_yaxis = []
     current_conc = concentrations[0]
     for i in range(100):
         concentration_line_xaxis.append(current_conc)
         concentration_lowerline_yaxis.append(lower_slope * current_conc + lower_intercept)
         concentration_higherline_yaxis.append(higher_slope * current_conc + higher_intercept)
     plt.scatter(concentration_line_xaxis, concentration_lowerline_yaxis, "r--")
     plt.scatter(concentration_line_xaxis, concentration_higherline_yaxis, "g--")
     cmc_x = []
     cmc_y = []
     cmc_x.append(cmc)
     cmc_y.append(cmc*lower_slope + lower_intercept)
     plt.scatter(cmc_x, cmc_y, "b")
     plt.show()


def test_main():
    cmc, lower_slope, lower_intercept, lower_rsqaured_optimized, higher_slope, higher_intercept, \
         higher_rsqaured_optimized, lower_linearregion_index_inclusive, \
              higher_linearregion_index_inclusive = determine_cmc(TEST_CONCENTRATIONS, TEST_CONDUCTIVITIES)
    plot_cmc(TEST_CONCENTRATIONS, TEST_CONDUCTIVITIES, cmc, lower_slope, lower_intercept, lower_rsqaured_optimized, higher_slope, \
     higher_intercept, higher_rsqaured_optimized, lower_linearregion_index_inclusive, higher_linearregion_index_inclusive)
    







        







