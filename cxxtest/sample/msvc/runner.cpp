/* Generated file, do not edit */

#ifndef CXXTEST_RUNNING
#define CXXTEST_RUNNING
#endif

#define _CXXTEST_HAVE_STD
#define _CXXTEST_HAVE_EH
#define _CXXTEST_ABORT_TEST_ON_FAIL
#include <cxxtest/TestListener.h>
#include <cxxtest/TestTracker.h>
#include <cxxtest/TestRunner.h>
#include <cxxtest/RealDescriptions.h>
#include <cxxtest/ParenPrinter.h>
#include <cxxtest/Win32Gui.h>

int main( int argc, char *argv[] ) {
 return CxxTest::GuiTuiRunner<CxxTest::Win32Gui, CxxTest::ParenPrinter>( argc, argv ).run();
}
#include "../gui/GreenYellowRed.h"

static GreenYellowRed suite_GreenYellowRed;

static CxxTest::List Tests_GreenYellowRed = { 0, 0 };
CxxTest::StaticSuiteDescription suiteDescription_GreenYellowRed( "../gui/GreenYellowRed.h", 11, "GreenYellowRed", suite_GreenYellowRed, Tests_GreenYellowRed );

static class TestDescription_GreenYellowRed_test_Start_green : public CxxTest::RealTestDescription {
public:
 TestDescription_GreenYellowRed_test_Start_green() : CxxTest::RealTestDescription( Tests_GreenYellowRed, suiteDescription_GreenYellowRed, 19, "test_Start_green" ) {}
 void runTest() { suite_GreenYellowRed.test_Start_green(); }
} testDescription_GreenYellowRed_test_Start_green;

static class TestDescription_GreenYellowRed_test_Green_again : public CxxTest::RealTestDescription {
public:
 TestDescription_GreenYellowRed_test_Green_again() : CxxTest::RealTestDescription( Tests_GreenYellowRed, suiteDescription_GreenYellowRed, 24, "test_Green_again" ) {}
 void runTest() { suite_GreenYellowRed.test_Green_again(); }
} testDescription_GreenYellowRed_test_Green_again;

static class TestDescription_GreenYellowRed_test_Now_yellow : public CxxTest::RealTestDescription {
public:
 TestDescription_GreenYellowRed_test_Now_yellow() : CxxTest::RealTestDescription( Tests_GreenYellowRed, suiteDescription_GreenYellowRed, 30, "test_Now_yellow" ) {}
 void runTest() { suite_GreenYellowRed.test_Now_yellow(); }
} testDescription_GreenYellowRed_test_Now_yellow;

static class TestDescription_GreenYellowRed_test_Cannot_go_back : public CxxTest::RealTestDescription {
public:
 TestDescription_GreenYellowRed_test_Cannot_go_back() : CxxTest::RealTestDescription( Tests_GreenYellowRed, suiteDescription_GreenYellowRed, 36, "test_Cannot_go_back" ) {}
 void runTest() { suite_GreenYellowRed.test_Cannot_go_back(); }
} testDescription_GreenYellowRed_test_Cannot_go_back;

static class TestDescription_GreenYellowRed_test_Finally_red : public CxxTest::RealTestDescription {
public:
 TestDescription_GreenYellowRed_test_Finally_red() : CxxTest::RealTestDescription( Tests_GreenYellowRed, suiteDescription_GreenYellowRed, 41, "test_Finally_red" ) {}
 void runTest() { suite_GreenYellowRed.test_Finally_red(); }
} testDescription_GreenYellowRed_test_Finally_red;

static class TestDescription_GreenYellowRed_test_Cannot_go_back_to_yellow : public CxxTest::RealTestDescription {
public:
 TestDescription_GreenYellowRed_test_Cannot_go_back_to_yellow() : CxxTest::RealTestDescription( Tests_GreenYellowRed, suiteDescription_GreenYellowRed, 47, "test_Cannot_go_back_to_yellow" ) {}
 void runTest() { suite_GreenYellowRed.test_Cannot_go_back_to_yellow(); }
} testDescription_GreenYellowRed_test_Cannot_go_back_to_yellow;

static class TestDescription_GreenYellowRed_test_Cannot_go_back_to_green : public CxxTest::RealTestDescription {
public:
 TestDescription_GreenYellowRed_test_Cannot_go_back_to_green() : CxxTest::RealTestDescription( Tests_GreenYellowRed, suiteDescription_GreenYellowRed, 53, "test_Cannot_go_back_to_green" ) {}
 void runTest() { suite_GreenYellowRed.test_Cannot_go_back_to_green(); }
} testDescription_GreenYellowRed_test_Cannot_go_back_to_green;

#include "../CreatedTest.h"

static CreatedTest *suite_CreatedTest = 0;

static CxxTest::List Tests_CreatedTest = { 0, 0 };
CxxTest::DynamicSuiteDescription<CreatedTest> suiteDescription_CreatedTest( "../CreatedTest.h", 14, "CreatedTest", Tests_CreatedTest, suite_CreatedTest, 21, 22 );

static class TestDescription_CreatedTest_test_nothing : public CxxTest::RealTestDescription {
public:
 TestDescription_CreatedTest_test_nothing() : CxxTest::RealTestDescription( Tests_CreatedTest, suiteDescription_CreatedTest, 24, "test_nothing" ) {}
 void runTest() { if ( suite_CreatedTest ) suite_CreatedTest->test_nothing(); }
} testDescription_CreatedTest_test_nothing;

#include "../DeltaTest.h"

static DeltaTest suite_DeltaTest;

static CxxTest::List Tests_DeltaTest = { 0, 0 };
CxxTest::StaticSuiteDescription suiteDescription_DeltaTest( "../DeltaTest.h", 7, "DeltaTest", suite_DeltaTest, Tests_DeltaTest );

static class TestDescription_DeltaTest_testSine : public CxxTest::RealTestDescription {
public:
 TestDescription_DeltaTest_testSine() : CxxTest::RealTestDescription( Tests_DeltaTest, suiteDescription_DeltaTest, 18, "testSine" ) {}
 void runTest() { suite_DeltaTest.testSine(); }
} testDescription_DeltaTest_testSine;

#include "../EnumTraits.h"

static EnumTraits suite_EnumTraits;

static CxxTest::List Tests_EnumTraits = { 0, 0 };
CxxTest::StaticSuiteDescription suiteDescription_EnumTraits( "../EnumTraits.h", 27, "EnumTraits", suite_EnumTraits, Tests_EnumTraits );

static class TestDescription_EnumTraits_test_Enum_traits : public CxxTest::RealTestDescription {
public:
 TestDescription_EnumTraits_test_Enum_traits() : CxxTest::RealTestDescription( Tests_EnumTraits, suiteDescription_EnumTraits, 30, "test_Enum_traits" ) {}
 void runTest() { suite_EnumTraits.test_Enum_traits(); }
} testDescription_EnumTraits_test_Enum_traits;

#include "../ExceptionTest.h"

static ExceptionTest suite_ExceptionTest;

static CxxTest::List Tests_ExceptionTest = { 0, 0 };
CxxTest::StaticSuiteDescription suiteDescription_ExceptionTest( "../ExceptionTest.h", 10, "ExceptionTest", suite_ExceptionTest, Tests_ExceptionTest );

static class TestDescription_ExceptionTest_testAssertion : public CxxTest::RealTestDescription {
public:
 TestDescription_ExceptionTest_testAssertion() : CxxTest::RealTestDescription( Tests_ExceptionTest, suiteDescription_ExceptionTest, 13, "testAssertion" ) {}
 void runTest() { suite_ExceptionTest.testAssertion(); }
} testDescription_ExceptionTest_testAssertion;

#include "../FixtureTest.h"

static FixtureTest suite_FixtureTest;

static CxxTest::List Tests_FixtureTest = { 0, 0 };
CxxTest::StaticSuiteDescription suiteDescription_FixtureTest( "../FixtureTest.h", 14, "FixtureTest", suite_FixtureTest, Tests_FixtureTest );

static class TestDescription_FixtureTest_test_strcpy : public CxxTest::RealTestDescription {
public:
 TestDescription_FixtureTest_test_strcpy() : CxxTest::RealTestDescription( Tests_FixtureTest, suiteDescription_FixtureTest, 28, "test_strcpy" ) {}
 void runTest() { suite_FixtureTest.test_strcpy(); }
} testDescription_FixtureTest_test_strcpy;

#include "../MessageTest.h"

static MessageTest suite_MessageTest;

static CxxTest::List Tests_MessageTest = { 0, 0 };
CxxTest::StaticSuiteDescription suiteDescription_MessageTest( "../MessageTest.h", 12, "MessageTest", suite_MessageTest, Tests_MessageTest );

static class TestDescription_MessageTest_testValues : public CxxTest::RealTestDescription {
public:
 TestDescription_MessageTest_testValues() : CxxTest::RealTestDescription( Tests_MessageTest, suiteDescription_MessageTest, 15, "testValues" ) {}
 void runTest() { suite_MessageTest.testValues(); }
} testDescription_MessageTest_testValues;

#include "../SimpleTest.h"

static SimpleTest suite_SimpleTest;

static CxxTest::List Tests_SimpleTest = { 0, 0 };
CxxTest::StaticSuiteDescription suiteDescription_SimpleTest( "../SimpleTest.h", 10, "SimpleTest", suite_SimpleTest, Tests_SimpleTest );

static class TestDescription_SimpleTest_testEquality : public CxxTest::RealTestDescription {
public:
 TestDescription_SimpleTest_testEquality() : CxxTest::RealTestDescription( Tests_SimpleTest, suiteDescription_SimpleTest, 13, "testEquality" ) {}
 void runTest() { suite_SimpleTest.testEquality(); }
} testDescription_SimpleTest_testEquality;

static class TestDescription_SimpleTest_testAddition : public CxxTest::RealTestDescription {
public:
 TestDescription_SimpleTest_testAddition() : CxxTest::RealTestDescription( Tests_SimpleTest, suiteDescription_SimpleTest, 21, "testAddition" ) {}
 void runTest() { suite_SimpleTest.testAddition(); }
} testDescription_SimpleTest_testAddition;

static class TestDescription_SimpleTest_TestMultiplication : public CxxTest::RealTestDescription {
public:
 TestDescription_SimpleTest_TestMultiplication() : CxxTest::RealTestDescription( Tests_SimpleTest, suiteDescription_SimpleTest, 27, "TestMultiplication" ) {}
 void runTest() { suite_SimpleTest.TestMultiplication(); }
} testDescription_SimpleTest_TestMultiplication;

static class TestDescription_SimpleTest_testComparison : public CxxTest::RealTestDescription {
public:
 TestDescription_SimpleTest_testComparison() : CxxTest::RealTestDescription( Tests_SimpleTest, suiteDescription_SimpleTest, 34, "testComparison" ) {}
 void runTest() { suite_SimpleTest.testComparison(); }
} testDescription_SimpleTest_testComparison;

static class TestDescription_SimpleTest_testTheWorldIsCrazy : public CxxTest::RealTestDescription {
public:
 TestDescription_SimpleTest_testTheWorldIsCrazy() : CxxTest::RealTestDescription( Tests_SimpleTest, suiteDescription_SimpleTest, 40, "testTheWorldIsCrazy" ) {}
 void runTest() { suite_SimpleTest.testTheWorldIsCrazy(); }
} testDescription_SimpleTest_testTheWorldIsCrazy;

static class TestDescription_SimpleTest_test_Failure : public CxxTest::RealTestDescription {
public:
 TestDescription_SimpleTest_test_Failure() : CxxTest::RealTestDescription( Tests_SimpleTest, suiteDescription_SimpleTest, 45, "test_Failure" ) {}
 void runTest() { suite_SimpleTest.test_Failure(); }
} testDescription_SimpleTest_test_Failure;

static class TestDescription_SimpleTest_test_TS_WARN_macro : public CxxTest::RealTestDescription {
public:
 TestDescription_SimpleTest_test_TS_WARN_macro() : CxxTest::RealTestDescription( Tests_SimpleTest, suiteDescription_SimpleTest, 51, "test_TS_WARN_macro" ) {}
 void runTest() { suite_SimpleTest.test_TS_WARN_macro(); }
} testDescription_SimpleTest_test_TS_WARN_macro;

#include "../TraitsTest.h"

static TestFunky suite_TestFunky;

static CxxTest::List Tests_TestFunky = { 0, 0 };
CxxTest::StaticSuiteDescription suiteDescription_TestFunky( "../TraitsTest.h", 53, "TestFunky", suite_TestFunky, Tests_TestFunky );

static class TestDescription_TestFunky_testPets : public CxxTest::RealTestDescription {
public:
 TestDescription_TestFunky_testPets() : CxxTest::RealTestDescription( Tests_TestFunky, suiteDescription_TestFunky, 56, "testPets" ) {}
 void runTest() { suite_TestFunky.testPets(); }
} testDescription_TestFunky_testPets;

#include <cxxtest/Root.cpp>
