#ifndef TRACE_TOGGLER_H
#define TRACE_TOGGLER_H

#include "trace-projections.h"
#include "converse.h"

// A converse msg to announce trace toggling
struct traceCmd
{
    char hdr[CmiMsgHeaderSizeBytes];
    bool isTracingNeeded;
};

// Fwd declaration
void traceCmdHandler(void *hdr);

/// A collection of routines that toggle projections tracing globally across the whole system
class traceToggler
{
    public:
        /// Converse msg handler that obeys an incoming trace cmd
        friend void traceCmdHandler(void *hdr)
        {
            traceCmd *cmd = (traceCmd*) hdr;
#ifdef LU_TRACING
            if (cmd->isTracingNeeded)
                traceBegin();
            else
                traceEnd();
#endif
            CmiFree(cmd);
        }

        /// initproc. Register a converse handler to accept the trace cmds
        static void registerHandler()
        {
#ifdef LU_TRACING
            traceCmdHandlerID = CmiRegisterHandler(traceCmdHandler);
#endif
        }

        /// Start projections tracing on this PE
        inline static void start()
        {
#ifdef LU_TRACING
            traceCmd *msg = (traceCmd*) CmiAlloc(sizeof(traceCmd));
            msg->isTracingNeeded = true;
            CmiSetHandler(msg, traceCmdHandlerID);
            CmiSyncBroadcastAllAndFree(sizeof(traceCmd), (char*)msg);
#endif
        }

        /// Stop projections tracing on this PE
        inline static void stop()
        {
#ifdef LU_TRACING
            traceCmd *msg = (traceCmd*) CmiAlloc(sizeof(traceCmd));
            msg->isTracingNeeded = false;
            CmiSetHandler(msg, traceCmdHandlerID);
            CmiSyncBroadcastAllAndFree(sizeof(traceCmd), (char*)msg);
#endif
        }

    private:
        static int traceCmdHandlerID;
};

#endif // TRACE_TOGGLER_H

