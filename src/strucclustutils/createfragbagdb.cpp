#include "Command.h"

#include "LocalParameters.h"
#include "DBReader.h"
#include "Debug.h"
#include "Util.h"
#include "FileUtil.h"
#include "Coordinate16.h"
#include "MultimerUtil.h"

#ifdef OPENMP
#include <omp.h>
#endif

int createfragbagdb(int argc, const char **argv, const Command &command) {
    // This function takes foldseek db as input, and creates a fragbag db based on the ca info
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, false, 0, MMseqsParameter::COMMAND_COMMON);

    return EXIT_SUCCESS;
}