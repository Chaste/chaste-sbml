#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace sbmltest
{
    /** Calculate the maximum of a vector of doubles.
     * @param vec The vector of doubles.
     * @return The maximum value.
     */
    inline double max(const std::vector<double> &vec);

    /** Calculate the mean of a vector of doubles.
     * @param vec The vector of doubles.
     * @return The mean value.
     */
    inline double mean(const std::vector<double> &vec);

    /** Calculate the minimum of a vector of doubles.
     * @param vec The vector of doubles.
     * @return The minimum value.
     */
    inline double min(const std::vector<double> &vec);

    /** Calculate the standard deviation of a vector of doubles.
     * @param vec The vector of doubles.
     * @return The standard deviation value.
     */
    inline double stdev(const std::vector<double> &vec);

    /** Calculate the qth quantile of a vector of doubles.
     * @param vec The vector of doubles, assumed to be sorted.
     * @param q The quantile to calculate (between 0 and 1).
     * @return The qth quantile value.
     */
    inline double quantile(const std::vector<double> &vec, double q);

    /** Calculate the variance of a vector of doubles.
     * @param vec The vector of doubles.
     * @return The variance value.
     */
    inline double variance(const std::vector<double> &vec);
} // namespace sbmltest

// max
inline double sbmltest::max(const std::vector<double> &vec)
{
    if (vec.empty())
    {
        throw std::invalid_argument("Cannot calculate maximum of an empty vector.");
    }
    return *std::max_element(vec.begin(), vec.end());
}

// mean
inline double sbmltest::mean(const std::vector<double> &vec)
{
    if (vec.empty())
    {
        throw std::invalid_argument("Cannot calculate mean of an empty vector.");
    }
    return std::accumulate(vec.begin(), vec.end(), 0.0) / vec.size();
}

// min
inline double sbmltest::min(const std::vector<double> &vec)
{
    if (vec.empty())
    {
        throw std::invalid_argument("Cannot calculate minimum of an empty vector.");
    }
    return *std::min_element(vec.begin(), vec.end());
}

// stdev
inline double sbmltest::stdev(const std::vector<double> &vec)
{
    return std::sqrt(sbmltest::variance(vec));
}

// quantile
inline double sbmltest::quantile(const std::vector<double> &vec, double q)
{
    if (vec.empty())
    {
        throw std::invalid_argument("Cannot calculate quantile of an empty vector.");
    }

    if (q < 0.0 || q > 1.0)
    {
        throw std::out_of_range("Quantile must be between 0.0 and 1.0");
    }

    std::vector<double> sorted_vec = vec;
    std::sort(sorted_vec.begin(), sorted_vec.end());

    size_t n = sorted_vec.size();
    size_t index = static_cast<size_t>(q * (n - 1));

    // Even number of elements
    if (n % 2 == 0)
    {
        return (sorted_vec[index] + sorted_vec[index - 1]) / 2.0;
    }

    // Odd number of elements
    return sorted_vec[index];
}

// variance
inline double sbmltest::variance(const std::vector<double> &vec)
{
    if (vec.size() < 2)
    {
        throw std::invalid_argument("Variance requires at least two data points.");
    }
    double mean_val = sbmltest::mean(vec);
    double accum = 0.0;
    for (double val : vec)
    {
        accum += (val - mean_val) * (val - mean_val);
    }
    return accum / (vec.size() - 1);
}
