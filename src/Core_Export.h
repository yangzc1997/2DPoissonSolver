// src/PoissonCore_Export.h
#ifndef POISSONCORE_EXPORT_H
#define POISSONCORE_EXPORT_H

#ifdef _WIN32
    #ifdef POISSONCORE_EXPORTS
        #define POISSONCORE_API __declspec(dllexport)
    #else
        #define POISSONCORE_API __declspec(dllimport)
    #endif
#else
    #define POISSONCORE_API
#endif

#endif // POISSONCORE_EXPORT_H