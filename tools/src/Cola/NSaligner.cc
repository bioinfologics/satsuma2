#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#include "NSaligner.h"


//=====================================================================
const AlignmentCola& NSaligner::align() {
  // No alignment can be performed if any of the query or target sequences are of 0 size
  if(!alignment.getTargetSeq().isize() || !alignment.getQuerySeq().isize()) { return alignment; } 
  double runtime = time(NULL);
  EditGraphDepth startPoint =  EditGraphDepth(editGraph.maxContigDepth, 0);
  double colaRuntimeFactor = traverseGraph(0, -1, editGraph.queryLen-1, editGraph.targetLen-1,
                          -1, startPoint);
  alignment.traceAlignment();
  runtime = time(NULL) - runtime;
  alignment.setRuntime(runtime);
  alignment.setRuntimeFactor(colaRuntimeFactor);
  return alignment;
}

double NSaligner::traverseGraph(int startRow, int startCol, int endRow, int endCol,
      int endDepth, const EditGraphDepth& prevCheckpointedCell) {
  double meanContigDepth = 0;
  EditGraphNode* currNode;
  // 1) The previously checkpointed cell should be set in its right place in the editGraph:
  // The place for this is given by the startRow and StartCol parameters
  editGraph.initCol(startCol, startRow, endRow);
  editGraph.getColumn(startCol)->setCell(startRow, prevCheckpointedCell); 

  // 2) Reset the bestScoredNode
  editGraph.updateBest(NULL);

  // 3) The middle column should be checkpointed
  int currCheckpointColIndex = (startCol + endCol)/2;

  // 4) Visit the graph nodes starting from the 
  // top most left node, moving vertically and to the right.
  // For this we start from the the column after the startCol
  // where the cell from previous calculations was set for restarting
  // calculations (hence the startCol+1 starting point below)
  for ( int col=startCol+1; col<=endCol; col++ ) {
    //Reset column
    editGraph.initCol(col, startRow, endRow); 
    //banded alignment - skip out-of-band cells
    int start = max(startRow, col-editGraph.bandWidth);
    int end   = min(endRow, col+editGraph.bandWidth);
    for ( int row=start; row<=end; row++ ) {
      int depth = 0;
      for ( depth; depth<=editGraph.maxContigDepth; depth++) {
        currNode = editGraph.getNode(row, col, depth);
        currNode->setCoords(row, col, depth); 
        visitNode(currNode, currCheckpointColIndex);
        if( currNode->getScore() == MINUS_INF && depth>2 ) { break; } //No need to search higher depths
       // (Step 4a) Update the current best node accordingly
        editGraph.updateBest(currNode);
        // (Step 4b) Set the best node for the current cell position
        if( currNode->getScore() > editGraph.getBestScoreAtRowCol(row, col)) { 
          editGraph.setBestNodeAtRowCol(row, col, currNode->getDepth()); 
        }
        // (Step 4c) For local alignment, if score is negative, set to zero 
        if(currNode->getScore()!=MINUS_INF && currNode->getScore()<0) {
          currNode->setScore(0);
        }
      }
      meanContigDepth += depth;
    }
    //Checkpoint middle column - keep for retrieving the checkpoint cell
    if(col == currCheckpointColIndex){ 
      editGraph.checkPoint(col, start, end);
    }
  }
  meanContigDepth /= (endRow*endCol);

  // 5) Check if reached end of recursion  
  if(startCol+1 == endCol) { 
    // If not at col 0, add to optimal path all the nodes on the column between the start and end row.
    if(endCol!=0) {
      for(int i=startRow+1; i<endRow; i++) { 
        alignment.addNodeToPath(editGraph.getBestNodeAtRowCol(i, endCol)); 
      }
    }
    return meanContigDepth;;
  }

  // 6) Track the ancestor node in the checkpointed column for the target node
  int optimalRow, optimalDepth;

  if(endDepth!=-1) { 
    optimalRow   = editGraph.getNode(endRow, endCol, endDepth)->getCPARow();
    optimalDepth = editGraph.getNode(endRow, endCol, endDepth)->getCPADepth();
  } else { //endDepth is set to -1 in the first run of the function as it's not known
    optimalRow   = editGraph.bestScoredNode.getCPARow();
    optimalDepth = editGraph.bestScoredNode.getCPADepth();
    //Change the coordinates of endNode so that the best node will become the end node
    endRow       = editGraph.bestScoredNode.getRow();
    endCol       = editGraph.bestScoredNode.getCol();
    endDepth     = editGraph.bestScoredNode.getDepth();
    // Add the targetNode to the path as this will not be added later
    alignment.addNodeToPath(&(editGraph.bestScoredNode));
    // If the best scoring node falls before the first checkpointed column:
    // Rerun function and do not continue to step 7-9 
    if(endCol<currCheckpointColIndex) {
      // Check to make sure endpoint hasn't reached start point (Very unlikely!)
      if ((startCol+1<endRow ) || (startCol+1==endRow && startRow<=endRow)) {
        traverseGraph(startRow, startCol, endRow, endCol,
           endDepth, prevCheckpointedCell);
      }   
      return meanContigDepth;
    }
  }
  // 7) First copy the checkpoint cell to be passed to next recursion
  EditGraphDepth currCheckpointCell = *(editGraph.checkpointCol.getCell(optimalRow));

  // 8) Then save node on the optimal path
  alignment.addNodeToPath(currCheckpointCell.getNode(optimalDepth));

  // 9) Continue recursion
  if ((startCol+1<currCheckpointColIndex ) || (startCol+1==currCheckpointColIndex && startRow<=optimalRow)) {
    traverseGraph(startRow, startCol, optimalRow, currCheckpointColIndex,
         optimalDepth, prevCheckpointedCell);
  }  
  if ((currCheckpointColIndex+1<endCol ) || (currCheckpointColIndex+1==endCol && optimalRow<=endRow)){
    traverseGraph(optimalRow, currCheckpointColIndex, endRow,
       endCol, endDepth, currCheckpointCell);
  } 

  return meanContigDepth;
}

void NSaligner::visitNode(EditGraphNode* currNode, int currCheckpointColIndex) {
  if (currNode->getDepth() == 0) {
    visitNodeContigZero(currNode, currCheckpointColIndex);
  } else {
    visitNodeContigNoneZero(currNode, currCheckpointColIndex);
  }
}

void NSaligner::visitNodeContigZero(EditGraphNode* currNode, int  currCheckpointColIndex) {
  // (Step 1) Find the score from all three possible neighbours (top,left,top-left)  
  int i = currNode->getRow();
  int j = currNode->getCol();
  int k = currNode->getDepth(); // Should be zero here
  double score1, score2, score3 = MINUS_INF;
  // Moving from the top or left neighbour, there is an indel, hence apply gap penalty
  score1 = editGraph.getBestScoreAtRowCol(i, j-1) + params.getGapOpenP();
  score2 = editGraph.getBestScoreAtRowCol(i-1, j) + params.getGapOpenP();  

  // Moving from the diagonal neighbour
  score3 = editGraph.getBestScoreAtRowCol(i-1, j-1) + params.getMismatchP();  
    
  // (Step 2) Find the maximum score and set the checkpoint ancestor coordinates
  if (( score3 >= score2 ) && ( score3 >= score1)) {
     currNode->setScore(score3);
    currNode->setCPACords(editGraph.getBestNodeAtRowCol(i-1, j-1),
      currCheckpointColIndex);
  } else if ( score2 > score1 ) {
    currNode->setScore(score2);
    currNode->setCPACords(editGraph.getBestNodeAtRowCol(i-1, j),
      currCheckpointColIndex);
  } else {
    currNode->setScore(score1);
    currNode->setCPACords(editGraph.getBestNodeAtRowCol(i, j-1),
      currCheckpointColIndex);
  }
} 

void NSaligner::visitNodeContigNoneZero(EditGraphNode* currNode, int  currCheckpointColIndex) {
  int k = currNode->getDepth();
  int i = currNode->getRow();
  int j = currNode->getCol();
  // In case there is no match:
  currNode->setScore(MINUS_INF);
  // Moving from the diagonal neighbour only if there's a match.
  if( getTargetSeq()[j] == getQuerySeq()[i] ) {
    double s = editGraph.getNode(i-1,j-1,k-1)->getScore();
    if(i*j==0 && k==1 && s==MINUS_INF) { s = 0; } // special case for first row/column
    if (s != MINUS_INF) { 
      //Scoring should become modular TODO 
//      currNode->setScore(s - exp(((k-1)-3)/1) + exp((k-3)/1));  
//      currNode->setScore(s - pow(1.5, ((k-1)-3)/1) + pow(1.5,(k-3)/1));  
      currNode->setScore(s - pow(((k-1)-0)/1.0, 3) + pow((k-0)/1.0, 3));  
      currNode->setCPACords(editGraph.getNode(i-1, j-1, k-1),
        currCheckpointColIndex);
    } 
  }
}

