#ifndef _EDITGRAPH_H_
#define _EDITGRAPH_H_

#include <limits>
#include <vector>
#include "../../analysis/DNAVector.h"

#define MINUS_INF  -numeric_limits<double>::max()

//==================================================================
/**
 * The nodes in the edit graph are represented by the following class.
 * Every node contains its coordinates (row, col), 
 *  a score for getting from it to the origin, a pointer
 * to the node which preceeded it in the route from origin to here,
 * and also the number of contiguous matches so far (Note that 
 * this means the number of matches immediately preceding the current 
 * node without having been interrupted by a mismatch or gap - this 
 * number can then be used to calculate the score of contiguous matches
 * and allocate a non-linear (exponential?) score to them.
 * Each node also contains the row index of the checkPointAncestor,
 * which is used for extracting the checkpoint row from the
 * best scored node once the give subgraph has been calculated.
 **/
class EditGraphNode
{
public:
  EditGraphNode():row(0), col(0), depth(0), score(MINUS_INF),
  CPArow(0), CPAdepth(0) {}
  int  getRow()               { return row; }
  int  getCol()               { return col; }
  int  getDepth()             { return depth; }
  double getScore()           { return score; }
  int  getCPARow()            { return CPArow; }
  int  getCPADepth()          { return CPAdepth; }
  void setScore(double s)     { score    = s; }
  void setRow(int r)          { row      = r; }
  void setCol(int c)          { col      = c; }
  void setDepth(int d)        { depth    = d; }
  void setCPARow(int r)       { CPArow   = r; }
  void setCPADepth(int d)     { CPAdepth = d; }

  /** Helper functions to set coordinates at the same time */
  void setCoords(int i, int j, int k) {
    setRow(i);
    setCol(j);
    setDepth(k);
  }
  /** This function sets the checkpoint ancestor coordinates to those passed via the parent node
   *  Note that this will be no-ops until the iteration has reached the midCol
   */
  void setCPACords(EditGraphNode* parent, int midCol);

private:
  int row;       /// The row index of the cell
  int col;       /// The column index of the cell
  int depth;     /// The number of contiguous matches
  double score;  /// The score for getting from the origin to this node
  int  CPArow;   /// Checkpoint Ancestor row
  int  CPAdepth; /// Checkpoint Ancestor depth 
};

//==================================================================
/** 
 * The third dimension of the edit graph are represented via this class.
 * Each instance contains all the nodes for a given row column at 
 * various contiguity depths from 0 to maxContigDepth. 
 * This class also holds the record for the node with the maximum score. 
 * This is then used in calculations instead of recalculating each time. 
 */
class EditGraphDepth
{
  public:
    EditGraphDepth(int maxContigDepth, int initScore=MINUS_INF): nodes(), bestNodeIndex(0) { 
      // Scores are all initialized to minus_inf but there are cases that requre otherwise
      getNode(0)->setScore(initScore);
    }
    ~EditGraphDepth() {}
    EditGraphNode* getNode(int k) { 
      if(getSize()<=k) {
        nodes.push_back(EditGraphNode());
      }
      return &nodes[k]; 
    } 
    EditGraphNode* getBestNode()  { return getNode(bestNodeIndex); }
    void setBestNode(int idx)     { bestNodeIndex = idx; } 
    double getBestScore() { 
      return getNode(bestNodeIndex)->getScore(); 
    }

    /** Return the number of nodes in the depth */
    int getSize() { return nodes.size(); }

    /** Used to initialise a cell depth for the next iteration */
    void init(int row, int col); 

  private:
    vector<EditGraphNode> nodes; /// Nodes with varying match contiguity depths.
    int bestNodeIndex;           /// Keep record to be used for scoring neighbouring cells
};

//==================================================================
/** 
 * Each column of the EditGraph contains a vector of EditGraphDepth cells 
 * This structure is used as the columns in the EditGraph class and also as
 * a container for the checkpoint columns. The column length corresponds
 * to the query length in the alignment + 1 (the addition is for the gap cell).
 */
class EditGraphColumn
{
public:
  EditGraphColumn(int qLen, int maxCD, int initScore=MINUS_INF):cells((qLen+1), EditGraphDepth(maxCD, initScore)) {} 

  ~EditGraphColumn() {}

  /**
   * Get a node at given coordinates.
   * @param[in]  The index in the query sequence
   * @return The requested node
   */
  EditGraphNode* getNode(int row, int depth) { return getCell(row)->getNode(depth); }

  /**
   * Get the EditGraphDepth at the given row position.
   * Note that row is incremented by 1 to cater for the buffer zone
   * @param[in]  The row position
   * @return The requested EditGraphDept
   */
  EditGraphDepth* getCell(int row) { return &(cells[row+1]); }

  /**
   * Used for places where a cell needs to be changed. E,g. At the beginning
   * of a recursion iteration where the cell is passed in from previous iteration.
   * Note that the +1 for the row index is to cater for the gap row.
   */
  void setCell(int row, const EditGraphDepth& cell) { cells[row+1] = cell; }
  void setCell(int row, EditGraphDepth* cell) { cells[row+1] = *cell; }

  /** Return the number of cells in the column */
  int getSize() { return cells.size(); }


private:
  vector< EditGraphDepth > cells; /// One dimension of the edit graph represented via a vector 
};


//==================================================================
// Forward Declaration
class NSaligner; 
class SWGAaligner; 
//==================================================================
/** 
 * The edit graph of the dynamic programming for the alignment
 * has two dimensions. Every cell contains a vector of nodes, 
 * each representing the node for the given cell position and
 * at a specific match contiguity depth. Using the check-point 
 * method only two columns of the graph are kept at any one time.
 * This is more space efficient than keeping the whole matrix
 * but it requires checkpointing (explained in algorithm methodology). 
 */
class EditGraph
{
  friend class NSaligner;
  friend class SWGAaligner;
public:
  EditGraph(int tLen, int qLen, int maxCD, int bandW):
    targetLen(tLen), queryLen(qLen), maxContigDepth(maxCD), 
    bandWidth(bandW), columns(2, EditGraphColumn(qLen+1, maxCD)),
    checkpointCol(qLen+1, maxCD), bestScoredNode() {
    //If bandwidth has not been provided, default is to run in unbanded mode
    if(bandWidth<0) { bandWidth = max(tLen, qLen); } 
  } 

  ~EditGraph() {}

  /**
   * Get a node at given coordinates.
   * @param[in]  The index in the target sequence
   * @param[in]  The index in the query sequence
   * @return The requested node
   */
  EditGraphNode* getNode(int row, int col, int depth) { 
    return getCell(row, col)->getNode(depth); }

  /**
   * Get the column that pertains to a column index
   * Note that the modular function is used based on the checkpointing algo 
   */
  EditGraphColumn* getColumn(int col) {
    return &(columns[((col+1)%2)]);  
  }

  /**
   * Get the best node from a set of nodes located at a specific cell position
   * @param[in] row: The row position of the cell
   * @param[in] col: The col position of the cell
   * @return The node with the best score at the given coordinates
   */
  EditGraphNode* getBestNodeAtRowCol(int row, int col) { return getCell(row, col)->getBestNode(); }
  void setBestNodeAtRowCol(int row, int col, int depthIdx) { 
    getCell(row, col)->setBestNode(depthIdx); 
  }
  double getBestScoreAtRowCol(int row, int col) { return getCell(row, col)->getBestScore(); }

  /**
   * Checks the current best and if the given node is better
   * the old one is deleted and replaced by the new one.
   * This function can also be used for resetting the bestNode
   * by passing Null as the bNode parameter
   */
  void updateBest(EditGraphNode* bNode); 
  
  /** 
   * Get the nodes at row,col 
   * Note that this function takes into consideration the fact
   * that only two columns of the EditGraph are being kept at
   * at any one time and hence it modulates the given column number by 2.
   * Also considers that the editGraph has been padded by an extra row
   * so passing indexes of -1 also works (all indexes are incremented by 1).
   */
  EditGraphDepth* getCell(int row, int col) { return getColumn(col)->getCell(row); }
 
  /**
   * Checks if a given cell (row, column) of the graph
   * is within the graphs bandWidth (for banded alignment)
   */
  bool isInBand(int row, int col) { return (abs(row-col)<=bandWidth); }
  bool isOnBandBorder(int row, int col) { return (abs(row-col)==(bandWidth+1)); }

protected:
  /** Used to initialize a column for the next iteration */
  void initCol(int col, int startRow, int endRow);

  /** Used to keep the checkpoint column updated so that the checkpointed node can be found. 
      This function does not copy the entire column as it might be that in the banded case and 
      also with smaller recursions the whole column is not needed.
  **/ 
  void checkPoint(int col, int maxStartRow, int minEndRow); 

  int targetLen;                   /// The number of characters ie target sequence
  int queryLen;                    /// The number of characters in the query sequence
  int maxContigDepth;              /// The maximum contiguity depth to be considered for scoring
  int bandWidth;                   /// The width of the band for banded smith-waterman
  vector<EditGraphColumn> columns; /// Two columns of the EditGraph kept at any one instance 
  EditGraphColumn checkpointCol;   /// Column used for checkpointing
  EditGraphNode bestScoredNode;    /// The node with the best score, used for tracing local alignment
};

#endif //_EDITGRAPH_H_
