//
// Created by stloufra on 11/19/24.
//
#pragma once

#include <fstream>
#include <iostream>
#include <string>
#include "../traits/LBMTraits.h"
#include <TNL/Timer.h>
#include <TNL/Logger.h>

#ifndef OUTPUTER_STATS
#define OUTPUTER_STATS

class outputerStats
{
    using RealType = LBMTraits::RealType;
    using DeviceType = LBMTraits::DeviceType;
    using LBMConstantsPointer = TNL::Pointers::SharedPointer<LBMConstants, DeviceType>;

public:
    template <typename SOL>
    static void Statistics(SOL& Solver, LBMConstantsPointer& Constants, bool withSys, bool withTur)
    {
        Logger logger(50, std::cout);
        if (withSys)
        {
            logger.writeSystemInformation(true);
        }
        logger.writeHeader("Timing of sections 1) Whole loop");
        Solver.timer_loop.writeLog(logger, 0);
        logger.writeSeparator();
        logger.writeHeader("Collision");
        Solver.timer_collision.writeLog(logger, 0);
        logger.writeSeparator();
        logger.writeHeader("Streaming");
        Solver.timer_streaming.writeLog(logger, 0);
        logger.writeSeparator();
        if (withTur)
        {
            logger.writeHeader("OmegaLes");
        }
        Solver.timer_LES.writeLog(logger, 0);
        logger.writeSeparator();
        logger.writeHeader("Bounce back");
        Solver.timer_bounceback.writeLog(logger, 0);
        logger.writeSeparator();
        logger.writeHeader("Moments Update");
        Solver.timer_momentsUpdate.writeLog(logger, 0);
        logger.writeSeparator();
        logger.writeHeader("Error Calculation");
        Solver.timer_err.writeLog(logger, 0);
        logger.writeSeparator();
        logger.writeHeader("Writting Output");
        Solver.timer_output.writeLog(logger, 0);
        logger.writeSeparator();
        logger.writeHeader("TimeDumping");
        Solver.timer_dumping.writeLog(logger, 0);
        logger.writeSeparator();
        logger.writeHeader("Time Averaging");
        Solver.timer_timeAvg.writeLog(logger, 0);

        auto MCells = (RealType)(Constants->dimX_int * Constants->dimY_int * Constants->dimZ_int) / 1e6;
        auto PerformTime = Solver.timer_collision.getRealTime()
            + Solver.timer_streaming.getRealTime()
            + Solver.timer_momentsUpdate.getRealTime()

            + Solver.timer_bounceback.getRealTime();

        if(withTur)
        {
            PerformTime += Solver.timer_LES.getRealTime();
        }

        auto iter = Constants->iterations;

        auto MLups = iter / PerformTime * MCells;

        printf("Performance - %f MLups \n", MLups);
    }
};

#endif