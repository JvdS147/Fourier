# Project: Fourier

CPP      = g++
CC       = gcc
OBJ = 3DCalculations.o AMS_Convert_flx2xyz.o AddClass.o AnalyseRings.o AnalyseTrajectory.o Angle.o AnisotropicDisplacementParameters.o Atom.o BackupOff0.o BagOfNumbers.o BasicMathsFunctions.o BondDetector.o CalculateBFDH.o ChebyshevBackground.o CheckFoundItem.o ChemicalFormula.o CollectionOfPoints.o Complex.o ConnectivityTable.o Constraints.o ConvexPolygon.o CoordinateFrame.o CopyTextFile.o CorrelationMatrix.o CrystalLattice.o CrystalStructure.o CyclicInteger.o Distance.o DoubleChecked.o DoubleWithESD.o DrunkardsWalk.o Eigenvalue.o Element.o EulerianAngles.o FileList.o FileName.o Finish_inp.o Fraction.o GenerateCombinations.o GeneratePowderCIF.o Histogram.o InpWriter.o LabelsAndShieldings.o MC_alkanes.o Mapping.o MathsFunctions.o Matrix3D.o MatrixFraction3D.o MillerIndices.o ModelBuilding.o MoleculeInCrystal.o NormalisedVector3D.o OneSudokuSlice.o OneSudokuSquare.o Plane.o PointGroup.o PowderMatchTable.o PowderPattern.o PowderPatternCalculator.o Pressure.o Quaternion.o RESULTSFile.o RandomNumberGenerator.o RandomQuaternionGenerator.o ReadCell.o ReadCif.o ReadCifOrCell.o ReadXSD.o ReadXYZ.o Refcode.o ReflectionList.o RunTests.o SetOfNumbers.o SimilarityAnalysis.o SkipBo.o SpaceGroup.o String2Fraction.o Sudoku.o SudokuSolver.o SymmetricMatrix3D.o SymmetryOperator.o TLS.o TLSWriter.o TOPAS.o Temperature.o Test3DCalculations.o TestAngle.o TestCalculateBFDH.o TestChebyshevBackground.o TestChemicalFormula.o TestComplex.o TestConstraints.o TestConvexPolygon.o TestCorrelationMatrix.o TestCrystalLattice.o TestCrystalStructure.o TestEulerianAngles.o TestFileName.o TestFraction.o TestGenerateCombinations.o TestMapping.o TestMaths.o TestMatrix3D.o TestMatrixFraction3D.o TestModelBuilding.o TestMoleculeInCrystal.o TestOneSudokuSlice.o TestOneSudokuSquare.o TestPowderMatchTable.o TestPowderPattern.o TestQuaternion.o TestRandomQuaternionGenerator.o TestReadCell.o TestReadXSD.o TestSetOfNumbers.o TestSort.o TestStack.o TestSudoku.o TestSudokuSolver.o TestSuite.o TestTLS.o TestTLS_ADPs.o TestTextFileReader_2.o TestTransSquareDependency.o TestTriangle.o TestTriangularPyramid.o TestUtilities.o TestVoidsFinder.o TextFileReader.o TextFileReader_2.o TextFileWriter.o TransSquareDependency.o Triangle.o TriangularPyramid.o Utilities.o Vector2D.o Vector3D.o Vector3DCalculations.o VoidsFinder.o Wavelength.o WriteCASTEPFile.o
LINKOBJ = $(OBJ)

CXXFLAGS = $(CXXINCS) -g -Ofast -Wfatal-errors #\
 #-DNORMALISE_HIGHEST_PEAK=${simulate_pattern_highest_peak} -DBACKGROUND_TOTAL_SIGNAL_NORMALISATION=${simulate_pattern_background} \
 #-DZERO_POINT_ERROR=${zero_point_error} -DFULL_WIDTH_HALF_MAX=${full_width_half_max} -DPREFERRED_ORIENTATION=${preferred_orientation}
CFLAGS   = $(INCS) -g -Ofast -Wfatal-errors
RM       = rm -f

all: Voids SimilarityMatrix UnitCellTransformation CentredToPrimitive

.PHONY: clean all clean-all

clean:
	$(RM) *.o

clean-all:
	$(RM) *.o Fourier Voids SimilarityMatrix

Fourier: Main.o $(OBJ)
	$(CPP) Main.o $(LINKOBJ) -o Fourier $(LIBS)

Voids: Voids.o $(OBJ)
	$(CPP) Voids.o $(LINKOBJ) -o Voids $(LIBS)

SimilarityMatrix: SimilarityMatrix.o $(OBJ)
	$(CPP) SimilarityMatrix.o $(LINKOBJ) -o SimilarityMatrix $(LIBS)

UnitCellTransformation: UnitCellTransformation.o $(OBJ)
	$(CPP) UnitCellTransformation.o $(LINKOBJ) -o UnitCellTransformation $(LIBS)

CentredToPrimitive: CentredToPrimitive.o $(OBJ)
	$(CPP) CentredToPrimitive.o $(LINKOBJ) -o CentredToPrimitive $(LIBS)

%.o: %.cpp
	$(CPP) -c $< -o $@ $(CXXFLAGS)
