///////////////////////// -*- C++ -*- /////////////////////////////
// ForwardJvtTool.cxx
// Implementation file for class ForwardJvtTool
// Author: Matt Klein<matthew.henry.klein@cern.ch>
///////////////////////////////////////////////////////////////////

// ForwardJvtTool includes
#include "ForwardJvtTool/ForwardJvtTool.h"

// Jet EDM
#include "xAODJet/JetAttributes.h"

// Shallow copy
#include "xAODCore/ShallowCopy.h"

    static SG::AuxElement::Decorator<char>  m_isHS("isJvtHS");
    static SG::AuxElement::Decorator<char>  m_isPU("isJvtPU");

  ///////////////////////////////////////////////////////////////////
  // Public methods:
  ///////////////////////////////////////////////////////////////////

  // Constructors
  ////////////////
  ForwardJvtTool::ForwardJvtTool(const std::string& name) :
    AsgTool(name)
  {
    declareProperty("OverlapDec",          m_orLabel            = ""                );
    declareProperty("OutputDec",           m_outLabel           = "passFJVT"        );
    declareProperty("EtaThresh",        m_etaThresh          = 2.5              );
    declareProperty("ForwardMinPt",        m_forwardMinPt          = 20e3              );
    declareProperty("ForwardMaxPt",        m_forwardMaxPt          = 50e3              );
    declareProperty("CentralMinPt",        m_centerMinPt          = 20e3              );
    declareProperty("CentralMaxPt",        m_centerMaxPt          = -1              );
    declareProperty("CentralJvtThresh",        m_centerJvtThresh          = 0.14              );
    declareProperty("CentralDrptThresh",       m_centerDrptThresh          = 0.26              );
    declareProperty("CentralMaxStochPt",          m_maxStochPt         = 35e3              );
    declareProperty("MaxDphi",          m_maxDphi         = 1.              );
  }

  // Destructor
  ///////////////
  ForwardJvtTool::~ForwardJvtTool()
  {}

  // Athena algtool's Hooks
  ////////////////////////////
  StatusCode ForwardJvtTool::initialize()
  {
    ATH_MSG_INFO ("Initializing " << name() << "...");
    if (m_orLabel!="")  Dec_OR = new SG::AuxElement::Decorator<char>(m_orLabel.Data());
    Dec_out = new SG::AuxElement::Decorator<char>(m_outLabel.Data());
    return StatusCode::SUCCESS;
  }

  StatusCode ForwardJvtTool::finalize()
  {
    ATH_MSG_INFO ("Finalizing " << name() << "...");
    return StatusCode::SUCCESS;
  }


  StatusCode ForwardJvtTool::processEvent(const xAOD::JetContainer* jets) {
    getPV();
    for(const auto& jetF : *jets) {
      (*Dec_out)(*jetF) = true;
      if (!forwardJet(jetF)) continue;
      for(const auto& jetC : *jets) {
        if (!centralJet(jetC)) continue;
        if (fabs(jetF->p4().DeltaPhi(-jetC->p4()))<m_maxDphi) (*Dec_out)(*jetF) = false;
      }
    }
    return StatusCode::SUCCESS;
  }

  bool ForwardJvtTool::passesFJVT(const xAOD::Jet *jet,const xAOD::JetContainer *jets) {
    if (!forwardJet(jet)) return true;
    getPV();
    for(const auto& jetC : *jets) {
      if (!centralJet(jetC)) continue;
      if (fabs(jet->p4().DeltaPhi(-jetC->p4()))<m_maxDphi) return false;
    }
    return true;
  }
  float ForwardJvtTool::getFJVT(const xAOD::Jet *jet,const xAOD::JetContainer *jets) {
    getPV();
    float bestmatch = 0;
    double tempdrpt = m_centerDrptThresh;
    m_centerDrptThresh = -1;
    for(const auto& jetC : *jets) {
      if (!centralJet(jetC)) continue;
      float current_drpt = jet->pt()>m_maxStochPt?0.99:getDrpt(jet);
      if (fabs(jet->p4().DeltaPhi(-jetC->p4()))<m_maxDphi && current_drpt>bestmatch) bestmatch = current_drpt;
    }
    m_centerDrptThresh = tempdrpt;
    return bestmatch;
  }

  bool ForwardJvtTool::forwardJet(const xAOD::Jet *jet) {
    if (fabs(jet->eta())<m_etaThresh) return false;
    if (jet->pt()<m_forwardMinPt || jet->pt()>m_forwardMaxPt) return false;
    return true;
  }

  bool ForwardJvtTool::centralJet(const xAOD::Jet *jet) {
    if (fabs(jet->eta())>m_etaThresh) return false;
    if (jet->pt()<m_centerMinPt || (m_centerMaxPt>0 && jet->pt()>m_centerMaxPt)) return false;
    if (Dec_OR && !(*Dec_OR)(*jet)) return false;
    float jvt;
    jet->getAttribute<float>("Jvt",jvt);
    if (jvt>m_centerJvtThresh) return false;
    if (jet->pt()<m_maxStochPt && getDrpt(jet)<m_centerDrptThresh) return false;
    return true;
  }

  float ForwardJvtTool::getDrpt(const xAOD::Jet *jet) {
    std::vector<float> sumpts;
    jet->getAttribute<std::vector<float> >("SumPtTrkPt500",sumpts);
    double firstVal = 0;
    double secondVal = 0;
    for (size_t i = 0; i < sumpts.size(); i++) {
      if (i==m_pvind) continue;
      if (sumpts[i]>firstVal) {
        secondVal = firstVal;
        firstVal = sumpts[i];
      } else if (sumpts[i]>secondVal) secondVal = sumpts[i];
    }
    return (firstVal-secondVal)/jet->pt();
  }

  void ForwardJvtTool::getPV() {
    const xAOD::VertexContainer *vxCont = 0;
    m_pvind = 0;
    if( evtStore()->retrieve(vxCont, "PrimaryVertices").isFailure() ) {
      ATH_MSG_WARNING("Unable to retrieve primary vertex container");
    } else if(vxCont->empty()) {
      ATH_MSG_WARNING("Event has no primary vertices!");
    } else {
      ATH_MSG_DEBUG("Successfully retrieved primary vertex container");
      for(const auto& vx : *vxCont) {
        if(vx->vertexType()==xAOD::VxType::PriVtx)
          {m_pvind = vx->index(); break;}
      }
    }
  }

  StatusCode ForwardJvtTool::tagTruth(const xAOD::JetContainer *jets,const xAOD::JetContainer *truthJets) {
    for(const auto& jet : *jets) {
      bool ishs = false;
      bool ispu = true;
      for(const auto& tjet : *truthJets) {
        if (tjet->p4().DeltaR(jet->p4())<0.3 && tjet->pt()>10e3) ishs = true;
        if (tjet->p4().DeltaR(jet->p4())<0.6) ispu = false;
      }
      m_isHS(*jet)=ishs;
      m_isPU(*jet)=ispu;
    }
    return StatusCode::SUCCESS;
  }

