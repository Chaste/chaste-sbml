
#include <fstream>

#include "OutputFileHandler.hpp"

#include "SbmlTestHelperFunctions.hpp"

void sbmltest::export_csv(const std::string& filename,
                          const std::vector<std::string>& var_names,
                          const OdeSolution& solutions)
{
    OutputFileHandler handler("");
    out_stream file = handler.OpenOutputFile(filename);

    // Write column headings
    (*file) << "time";
    for (unsigned i = 0; i < var_names.size(); i++)
    {
        (*file) << ", " << var_names[i];
    }
    (*file) << "\n"
            << std::flush;

    // Write data
    for (unsigned i = 0; i < solutions.rGetTimes().size(); i++)
    {
        (*file) << solutions.rGetTimes()[i];
        for (unsigned j = 0; j < solutions.rGetSolutions()[i].size(); j++)
        {
            (*file) << ", " << solutions.rGetSolutions()[i][j];
        }
        (*file) << "\n"
                << std::flush;
    }
    file->close();
}
