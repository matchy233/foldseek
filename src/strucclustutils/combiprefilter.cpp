#include "Command.h"
#include "Debug.h"

extern int prefilter(int argc, const char **argv, const Command &command);
int combiprefilter(int argc, const char **argv, const Command &command)
{
    Debug(Debug::INFO) << "Running combined prefilter instead of old prefilter\n";
    return prefilter(argc, argv, command);
}