///////////////////////// -*- C++ -*- /////////////////////////////
// ForwardJvtTool.h
// Header file for class ForwardJvtTool
// Author: Matt Klein<matthew.henry.klein@cern.ch>
///////////////////////////////////////////////////////////////////
#ifndef FORWARDJVTTOOL_JVT_FORWARDJVTTOOL_H
#define FORWARDJVTTOOL_JVT_FORWARDJVTTOOL_H 1

// STL includes
#include <string>

// FrameWork includes
#include "AsgTools/ToolHandle.h"
#include "AsgTools/AsgTool.h"

// EDM includes
#include "xAODJet/JetContainer.h"

  class ForwardJvtTool
  : virtual public asg::AsgTool

  {
    // This macro defines the constructor with the interface declaration
    //ASG_TOOL_CLASS(ForwardJvtTool, IForwardJvtTool)

    ///////////////////////////////////////////////////////////////////
    // Public methods:
    ///////////////////////////////////////////////////////////////////
  public:

    // Copy constructor:

    /// Constructor with parameters:
    ForwardJvtTool(const std::string& name);

    /// Destructor:
    virtual ~ForwardJvtTool();

    // Athena algtool's Hooks
    StatusCode  initialize();
    StatusCode  finalize();

    StatusCode processEvent(const xAOD::JetContainer* jets);
    bool passesFJVT(const xAOD::Jet *jet,const xAOD::JetContainer *jets);
    float getFJVT(const xAOD::Jet *jet,const xAOD::JetContainer *jets);
    bool forwardJet(const xAOD::Jet *jet);
    bool centralJet(const xAOD::Jet *jet);
    float getDrpt(const xAOD::Jet *jet);

    static StatusCode tagTruth(const xAOD::JetContainer *jets,const xAOD::JetContainer *truthJets);

  private:

    TString m_orLabel;
    TString m_outLabel;
    double m_etaThresh;
    double m_forwardMinPt;
    double m_forwardMaxPt;
    double m_centerMinPt;
    double m_centerMaxPt;
    double m_centerJvtThresh;
    double m_centerDrptThresh;
    double m_maxStochPt;
    double m_maxDphi;
    size_t m_pvind;
    SG::AuxElement::Decorator<char>* Dec_OR = NULL;
    SG::AuxElement::Decorator<char>* Dec_out = NULL;
    void getPV();

    /// Default constructor:
    ForwardJvtTool();

  };
#endif //> !FORWARDJVTTOOL_JVT_FORWARDJVTTOOL_H
