#ifndef _ALIGNER_PARAMS_H_
#define _ALIGNER_PARAMS_H_

//=====================================================================

/**
 * Enumeration specifying the different aligner types available within Cola
 * NSGAaligner (specified by NSGA:1), 
 * NSaligner (specified by NS:2)
 * SWGAaligner (specified by SWGA:3)
 * SWaligner (specified by SW:4)
 */
enum AlignerType { UNUSED, NSGA, NS, SWGA, SW };

/**
 * Object encapsulating necessary parameters and type identifier for setting up an aligner
 */
class AlignerParams {
public:
  // Default Ctor
  AlignerParams():bandWidth(-1), alignerType(NSGA), useAlignerDef(true),
      gapOpenP(0), mismatchP(0), gapExtP(0) { setDefaults(); }
  // Ctor 2
  AlignerParams(int bandW):bandWidth(bandW), alignerType(NSGA), useAlignerDef(true),
      gapOpenP(0), mismatchP(0), gapExtP(0) { setDefaults(); }
  // Ctor 2
  AlignerParams(int bandW, AlignerType type):bandWidth(bandW), alignerType(type), useAlignerDef(true),
      gapOpenP(0), mismatchP(0), gapExtP(0) { setDefaults(); }
  // Ctor 3
  AlignerParams(int bandW, AlignerType type, int goPen, int mmPen,
       int gePen):bandWidth(bandW), alignerType(type), useAlignerDef(false),
        gapOpenP(goPen), mismatchP(mmPen), gapExtP(gePen) {}

// Setters
  void setType(AlignerType at)   { alignerType = at;  }
  void setGapOpenP(int gop)      { gapOpenP    = gop; }
  void setMismatchP(int mp)      { mismatchP   = mp;  }
  void setGapExtP(int gep)       { gapExtP     = gep; }
  void setBandWidth(int bw)      { bandWidth   = bw;  }

// Getters
  AlignerType getType()const   { return alignerType; }
  int  getGapOpenP()const      { return gapOpenP; }
  int  getMismatchP()const     { return mismatchP; }
  int  getGapExtP()const       { return gapExtP; }
  bool useDefaults()const      { return useAlignerDef; }
  int  getBandWidth()const     { return bandWidth; }

private:
  /** Set the param defaults if they haven't been given in the constructor */
  void setDefaults() {
    switch(alignerType) {
      case NSGA:
        gapOpenP  = -200;
	mismatchP = -8;
        gapExtP   = -20;
	break;
      case NS:
        gapOpenP  = -100;
        mismatchP = -8;
        gapExtP   = -100; // Same as gapOpen
	break;
      case SWGA:
        gapOpenP  = -6;
        mismatchP = -1;
        gapExtP   = -1; 
	break;
      case SW:
        gapOpenP  = -1;
        mismatchP = -1;
        gapExtP   = -1; // Same as gapOpen
	break;
      default:
        gapOpenP  = -200;
	mismatchP = -8;
        gapExtP   = -20;
  }	
}

  int bandWidth;           /// The bandwidth for banded alignment. Use -1 for unbanded alignment
  AlignerType alignerType; /// Type of the underlying aligner to use for refinement
  bool useAlignerDef;      /// Flag specifying whether to use aligner defaults or not
  int gapOpenP;            /// Gap Open Penalty
  int mismatchP;           /// Mismatch penalty
  int gapExtP;             /// Gap extension penalty
};

#endif //_ALIGNER_PARAMS_H_
