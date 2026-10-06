
#pragma once

// Declarations for logging helpers provided by the main module.
// Implementations are in src/ServerMain.cpp in this Linux port.

void ErrorToLog(const char *message);
void DebugToLog(const char *message);
void EitherToLog(const char *prefix, const char *message, bool saveErr = false);
