/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package macro;

import java.util.*;
import java.io.*;
import java.nio.*;
import star.flow.*;
import star.base.neo.*;   // for use with version >= 1.09
import star.common.*;
import star.coremodule.SimulationProcessObject;
import star.coremodule.services.*;
import star.saturb.*;
import star.keturb.*;
import star.kwturb.*;
//import star.segregatedenergy.*;
//import star.segregatedflow.*;
import star.coupledflow.*;

/**
 *
 * @author tjp644
 */
public class ClearData {
    Simulation sim = 
      getActiveSimulation();
    String simName = sim.getPresentationName();
    public void execute() {
        sim.clearSolution();
        sim.getMeshPipelineController().clearGeneratedMeshes();
        sim.saveState(simName+"empty.sim");
    }
}
