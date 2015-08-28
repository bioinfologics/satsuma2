export CPP="icpc -g -std=c++0x -O3 -fno-math-errno"
$CPP -o HomologyByXCorrSlave analysis/HomologyByXCorrSlave.cc analysis/ProbTable.cc analysis/DNAVector.cc base/FileParser.cc util/mutil.cc analysis/CrossCorr.cc analysis/CodonTranslate.cc analysis/SeqChunk.cc base/StringUtil.cc analysis/AlignProbability.cc -I . -lpthread &
$CPP -o SatsumaSynteny2 analysis/SatsumaSynteny2.cc -I . analysis/DNAVector.cc base/FileParser.cc util/mutil.cc analysis/CodonTranslate.cc analysis/SeqChunk.cc base/StringUtil.cc analysis/AlignProbability.cc analysis/WorkQueue.cc analysis/MatchDynProg.cc util/SysTime.cc analysis/SequenceMatch.cc analysis/GridSearch.cc -lpthread &
$CPP -o MergeXCorrMatches analysis/MergeXCorrMatches.cc base/StringUtil.cc analysis/DNAVector.cc base/FileParser.cc analysis/XCorrDynProg.cc analysis/SequenceMatch.cc util/mutil.cc  -I . &
$CPP -o MatchesByFeature tools/MatchesByFeature.cc &
$CPP -o KMatch kmatch/KMatch.cc -lpthread -I . &
$CPP -o matchdump kmatch/matchdump.cc -I . &
$CPP -o ChainMatches analysis/ChainMatches.cc base/StringUtil.cc analysis/DNAVector.cc base/FileParser.cc analysis/XCorrDynProg.cc analysis/MatchDynProg.cc analysis/SequenceMatch.cc util/mutil.cc  -I . &

wait
