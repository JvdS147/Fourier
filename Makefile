# Project: Fourier

CPP      = g++
CC       = gcc
OBJ = Main.o 3DCalculations.o CrystalLattice.o MathFunctions.o Refcode.o TestCorrelationMatrix.o TestSuite.o AMS_Convert_flx2xyz.o CrystalStructure.o Matrix3D.o ReflectionList.o TestCrystalLattice.o TestTLS_ADPs.o AddClass.o CyclicInteger.o MillerIndices.o RunTests.o TestCrystalStructure.o TestTransSquareDependency.o AnalyseRings.o Distance.o ModelBuilding.o SetOfNumbers.o TestFileName.o TestTriangle.o AnalyseTrajectory.o DoubleChecked.o MoleculeInCrystal.o SimilarityAnalysis.o TestFraction.o TestTriangularPyramid.o Angle.o DoubleWithESD.o NormalisedVector3D.o SkipBo.o TestGenerateCombinations.o TestUtilities.o AnisotropicDisplacementParameters.o DrunkardsWalk.o OneSudokuSlice.o SpaceGroup.o TestMatrix3D.o TestVoidsFinder.o Atom.o Eigenvalue.o OneSudokuSquare.o String2Fraction.o TestModelBuilding.o TextFileReader.o BackupOff0.o Element.o Plane.o Sudoku.o TestMoleculeInCrystal.o TextFileReader_2.o BagOfNumbers.o FileList.o PointGroup.o SudokuSolver.o TestOneSudokuSlice.o TextFileWriter.o BondDetector.o FileName.o PowderMatchTable.o SymmetricMatrix3D.o TestOneSudokuSquare.o TransSquareDependency.o CalculateBFDH.o Finish_inp.o PowderPattern.o SymmetryOperator.o TestPowderMatchTable.o Triangle.o ChebyshevBackground.o Fraction.o PowderPatternCalculator.o TLSWriter.o TestQuaternion.o TriangularPyramid.o CheckFoundItem.o GenerateCombinations.o Pressure.o TOPAS.o TestRandomQuaternionGenerator.o Utilities.o ChemicalFormula.o GeneratePowderCIF.o Quaternion.o Temperature.o TestReadXSD.o Vector3D.o CollectionOfPoints.o Histogram.o RandomNumberGenerator.o Test3DCalculations.o TestSetOfNumbers.o Vector3DCalculations.o ConnectivityTable.o InpWriter.o RandomQuaternionGenerator.o TestAngle.o TestSort.o VoidsFinder.o ConvexPolygon.o LabelsAndShieldings.o ReadCif.o TestCalculateBFDH.o TestStack.o Wavelength.o CopyTextFile.o ReadXSD.o TestChebyshevBackground.o TestSudoku.o WriteCASTEPFile.o CorrelationMatrix.o MathConstants.o ReadXYZ.o TestConvexPolygon.o TestSudokuSolver.o
LINKOBJ = Main.o 3DCalculations.o CrystalLattice.o MathFunctions.o Refcode.o TestCorrelationMatrix.o TestSuite.o AMS_Convert_flx2xyz.o CrystalStructure.o Matrix3D.o ReflectionList.o TestCrystalLattice.o TestTLS_ADPs.o AddClass.o CyclicInteger.o MillerIndices.o RunTests.o TestCrystalStructure.o TestTransSquareDependency.o AnalyseRings.o Distance.o ModelBuilding.o SetOfNumbers.o TestFileName.o TestTriangle.o AnalyseTrajectory.o DoubleChecked.o MoleculeInCrystal.o SimilarityAnalysis.o TestFraction.o TestTriangularPyramid.o Angle.o DoubleWithESD.o NormalisedVector3D.o SkipBo.o TestGenerateCombinations.o TestUtilities.o AnisotropicDisplacementParameters.o DrunkardsWalk.o OneSudokuSlice.o SpaceGroup.o TestMatrix3D.o TestVoidsFinder.o Atom.o Eigenvalue.o OneSudokuSquare.o String2Fraction.o TestModelBuilding.o TextFileReader.o BackupOff0.o Element.o Plane.o Sudoku.o TestMoleculeInCrystal.o TextFileReader_2.o BagOfNumbers.o FileList.o PointGroup.o SudokuSolver.o TestOneSudokuSlice.o TextFileWriter.o BondDetector.o FileName.o PowderMatchTable.o SymmetricMatrix3D.o TestOneSudokuSquare.o TransSquareDependency.o CalculateBFDH.o Finish_inp.o PowderPattern.o SymmetryOperator.o TestPowderMatchTable.o Triangle.o ChebyshevBackground.o Fraction.o PowderPatternCalculator.o TLSWriter.o TestQuaternion.o TriangularPyramid.o CheckFoundItem.o GenerateCombinations.o Pressure.o TOPAS.o TestRandomQuaternionGenerator.o Utilities.o ChemicalFormula.o GeneratePowderCIF.o Quaternion.o Temperature.o TestReadXSD.o Vector3D.o CollectionOfPoints.o Histogram.o RandomNumberGenerator.o Test3DCalculations.o TestSetOfNumbers.o Vector3DCalculations.o ConnectivityTable.o InpWriter.o RandomQuaternionGenerator.o TestAngle.o TestSort.o VoidsFinder.o ConvexPolygon.o LabelsAndShieldings.o ReadCif.o TestCalculateBFDH.o TestStack.o Wavelength.o CopyTextFile.o ReadXSD.o TestChebyshevBackground.o TestSudoku.o WriteCASTEPFile.o CorrelationMatrix.o MathConstants.o ReadXYZ.o TestConvexPolygon.o TestSudokuSolver.o

BIN      = Fourier
CXXFLAGS = $(CXXINCS) -Ofast -Wfatal-errors
CFLAGS   = $(INCS) -Ofast -Wfatal-errors
RM       = rm -f

all: $(BIN)

.PHONY: clean all

clean:
	$(RM) $(OBJ) $(BIN)

$(BIN): $(OBJ)
	$(CPP) $(LINKOBJ) -o $(BIN) $(LIBS)

$(OBJ): %.o: %.cpp
	$(CPP) -c $< -o $@ $(CXXFLAGS)
