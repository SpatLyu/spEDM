#ifndef DeLongPlacements_H
#define DeLongPlacements_H

#include <vector>
#include <algorithm>
#include <utility>
#include <string>
#include "spEDMDataStruct.h"

/**
 * @brief Computes DeLong placements for ROC analysis.
 *
 * This function implements the DeLong method to calculate placements for cases and controls
 * in the context of ROC (Receiver Operating Characteristic) curve analysis. The method is
 * used for non-parametric estimation of the area under the curve (AUC) and supports handling
 * tied values efficiently.
 *
 * The function accepts two sets of data points: `cases` and `controls`. Depending on the
 * specified `direction` (either ">" or "<"), it processes the data accordingly by inverting
 * the sign of the input values if necessary.
 *
 * @param cases A vector of numeric values representing positive cases.
 * @param controls A vector of numeric values representing negative controls.
 * @param direction A string indicating the comparison direction. If set to ">", values are inverted.
 * @return DeLongPlacementsRes A structure containing:
 *   - `theta`: The estimated AUC value.
 *   - `X`: A vector of normalized placement values for cases.
 *   - `Y`: A vector of normalized placement values for controls.
 */
DeLongPlacementsRes CppDeLongPlacements(const std::vector<double>& cases,
                                        const std::vector<double>& controls,
                                        const std::string& direction);

#endif // DeLongPlacements_H
