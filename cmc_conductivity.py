import math
import statistics
import matplotlib.pyplot as plt
import csv

def determine_cmc(concentrations, conductivities):
    lower_rsqaured_optimized = 0
    lower_linearregion_index_inclusive = 2
    for i in range(3, len(concentrations)):
        lower_rsquared = statistics.correlation(concentrations[0:i], conductivities[0:i]) ** 2
        if lower_rsquared >= lower_rsqaured_optimized:
            lower_rsqaured_optimized = lower_rsquared
            lower_linearregion_index_inclusive = i - 1
    higher_rsqaured_optimized = 0
    higher_linearregion_index_inclusive = lower_linearregion_index_inclusive + 1
    for i in range(3, len(concentrations)):
        higher_rsquared = statistics.correlation(concentrations[len(concentrations)-(i+1):len(concentrations)-1], 
        conductivities[len(conductivities)-(i+1):len(conductivities)-1]) ** 2
        if higher_rsquared >=  higher_rsqaured_optimized:
            higher_rsqaured_optimized = higher_rsquared
            higher_linearregion_index_inclusive = len(concentrations) - i
    lower_slope, lower_intercept = statistics.linear_regression(concentrations[0:lower_linearregion_index_inclusive+1], \
         conductivities[0:lower_linearregion_index_inclusive+1])
    higher_slope, higher_intercept = statistics.linear_regression(concentrations[higher_linearregion_index_inclusive:len(concentrations)], \
         conductivities[higher_linearregion_index_inclusive:len(conductivities)])
    cmc = (higher_intercept - lower_intercept) / (lower_slope - higher_slope)
    return cmc, lower_slope, lower_intercept, lower_rsqaured_optimized, higher_slope, higher_intercept, higher_rsqaured_optimized, lower_linearregion_index_inclusive, higher_linearregion_index_inclusive

def determine_dgmicellization(temp, num_chargedgroups, charge_pergroup, num_tails, charge_percounterion, lower_slope, higher_slope, cmc):
    beta = (lower_slope-higher_slope)/lower_slope
    term1 = 8.314 * temp * (1/num_tails + beta * (num_chargedgroups/num_tails) \
    * math.fabs(charge_pergroup/charge_percounterion)) * math.log(cmc)
    term2 = 8.314 * temp * (beta * ((num_chargedgroups/num_tails)*math.fabs(charge_pergroup/charge_percounterion)) \
        * math.log(((num_chargedgroups/num_tails)*math.fabs(charge_pergroup/charge_percounterion))) \
        - math.log(num_tails)/num_tails)
    return term1 + term2

def plot_cmc(concentrations, conductivities, cmc, lower_slope, lower_intercept, lower_rsqaured_optimized, higher_slope, \
     higher_intercept, higher_rsqaured_optimized, lower_linearregion_index_inclusive, higher_linearregion_index_inclusive, notes):
     plt.scatter(concentrations,conductivities,color="k")
     plt.scatter(concentrations[0:lower_linearregion_index_inclusive+1], conductivities[0:lower_linearregion_index_inclusive+1], color="r")
     plt.scatter(concentrations[higher_linearregion_index_inclusive:len(concentrations)], conductivities[higher_linearregion_index_inclusive:len(concentrations)], color="g")
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
         current_conc += concentration_stepsize
     plt.plot(concentration_line_xaxis, concentration_lowerline_yaxis, color="r", linestyle="dashed")
     plt.plot(concentration_line_xaxis, concentration_higherline_yaxis, color="g", linestyle="dashed")
     cmc_x = []
     cmc_y = []
     cmc_x.append(cmc)
     cmc_y.append(cmc*lower_slope + lower_intercept)
     plt.scatter(cmc_x, cmc_y, color="b")
     plt.xlabel("Surfactant Concentration (mM)")
     plt.ylabel("Conductivity (" + u"\u03bcS/cm)")
     plt.title(notes)
     plt.show()


def filter_measurements(concentrations_init, conductivities_init):
    bad_indexes = []
    for i in range(len(concentrations_init)):
        try:
            test_conc = float(concentrations_init[i])
            test_cond = int(conductivities_init[i])
        except:
            bad_indexes.append(i)
    concentrations = []
    conductivities = []
    for i in range(len(concentrations_init)):
        if i not in bad_indexes:
            concentrations.append(concentrations_init[i])
            conductivities.append(conductivities_init[i])
    return concentrations, conductivities

def read_csvdata(csv_filename, mode_readsimple=True):
    if mode_readsimple:
        concentrations = []
        conductivities = []
        with open (csv_filename + ".csv", mode="r") as input_csvfile:
            input_csvreader = csv.reader(input_csvfile)
            row_counter = 0
            for row in input_csvreader:
                if row_counter != 0:
                    try:
                        concentrations.append(float(row[0]))
                        conductivities.append(float(row[1]))
                    except:
                        pass
                row_counter += 1
        input_csvfile.close()
    return concentrations, conductivities


def write_csvdata(csv_filename, concentrations, conductivities):
    with open (csv_filename + ".csv", mode="w") as output_csvfile:
        output_csvwriter = csv.writer(output_csvfile)
        output_csvwriter.writerow(["Surfactant Concentration (mM)", "Conductivity (uS/cm)"])
        for i in range(len(concentrations)):
            output_csvwriter.writerow([str(concentrations[i]),str(conductivities[i])])
    output_csvfile.close()
    return

def write_txtdata(txt_filename, concentrations, conductivities, notes, cmc, lower_slope, lower_intercept, \
     lower_rsqaured_optimized, higher_slope, higher_intercept, higher_rsqaured_optimized, lower_linearregion_index_inclusive, \
         higher_linearregion_index_inclusive):
    with open(txt_filename + ".txt", "w") as output_txtfile:
        output_txtfile.write("NOTES FOR RUN: \n" + notes + "\n\n")
        output_txtfile.write("SUMMARY OF DATA ANALYSIS FOR RUN: \n")
        output_txtfile.write("CMC: " + str(cmc) + " mM\n")
        output_txtfile.write("Lower Line: y = " + str(lower_slope) + "x + " + str(lower_intercept) + " (R^2 = " + str(lower_rsqaured_optimized) + ")\n")
        output_txtfile.write("Num. Points Used: " + str(lower_linearregion_index_inclusive + 1) + "\n")
        output_txtfile.write("Higher Line: y = " + str(higher_slope) + "x + " + str(higher_intercept) + " (R^2 = " + str(higher_rsqaured_optimized) + ")\n")
        output_txtfile.write("Num. Points Used: " + str(len(concentrations) - higher_linearregion_index_inclusive) + "\n\n\n")
        output_txtfile.write("DATA FOR RUN: \n")
        output_txtfile.write("Concentrations (mM), Conductivity (uS/cm)\n")
        for i in range(len(concentrations)):
            output_txtfile.write(str(concentrations[i]) + ", " + str(conductivities[i]) + "\n")
    output_txtfile.close()




def test_main(notes, input_filename, output_filename, temp, num_chargedgroups, charge_pergroup, num_tails, charge_percounterion):
    TEST_CONCENTRATIONS, TEST_CONDUCTIVITIES = read_csvdata(input_filename)
    concentrations, conductivities = filter_measurements(TEST_CONCENTRATIONS, TEST_CONDUCTIVITIES)
    cmc, lower_slope, lower_intercept, lower_rsqaured_optimized, higher_slope, higher_intercept, \
         higher_rsqaured_optimized, lower_linearregion_index_inclusive, \
              higher_linearregion_index_inclusive = determine_cmc(concentrations, conductivities)
    dg_micellization = determine_dgmicellization(temp, num_chargedgroups, charge_pergroup, num_tails, charge_percounterion, lower_slope, higher_slope, cmc)
    print(dg_micellization)
    write_csvdata(output_filename, concentrations, conductivities)
    write_txtdata(output_filename, concentrations, conductivities, notes, cmc, lower_slope, lower_intercept, \
     lower_rsqaured_optimized, higher_slope, higher_intercept, higher_rsqaured_optimized, lower_linearregion_index_inclusive, \
         higher_linearregion_index_inclusive)
    plot_cmc(concentrations, conductivities, cmc, lower_slope, lower_intercept, lower_rsqaured_optimized, higher_slope, \
     higher_intercept, higher_rsqaured_optimized, lower_linearregion_index_inclusive, higher_linearregion_index_inclusive, notes)


if __name__ == "__main__":
    test_main("Undecyl Alanine in the Presence of Ethylenediamine Counterion\npH 8 at 298 K", "input", "output", 298, 1, -1, 1, 1)
    







        







