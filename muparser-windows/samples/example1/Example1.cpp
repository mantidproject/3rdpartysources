//---------------------------------------------------------------------------
//
//                 __________                                      
//    _____   __ __\______   \_____  _______  ______  ____ _______ 
//   /     \ |  |  \|     ___/\__  \ \_  __ \/  ___/_/ __ \\_  __ \ 
//  |  Y Y  \|  |  /|    |     / __ \_|  | \/\___ \ \  ___/ |  | \/
//  |__|_|  /|____/ |____|    (____  /|__|  /____  > \___  >|__|   
//        \/                       \/            \/      \/        
//  (C) 2004-2007 Ingo Berg
//
//  Example 1 - using the parser as a static library
//
//---------------------------------------------------------------------------

#include "muParserTest.h"

/** \brief This macro will enable mathematical constants like M_PI. */
#define _USE_MATH_DEFINES		

#include <cstdlib>
#include <cstring>
#include <cmath>
#include <string>
#include <iostream>

#include "muParser.h"

#if defined( USINGDLL ) && defined( _WIN32 )
#error This sample can be used only with STATIC builds of muParser (on win32)
#endif


using namespace std;
using namespace mu;

// Operator callback functions
value_type Mega(value_type a_fVal) { return a_fVal * 1e6; }
value_type Milli(value_type a_fVal) { return a_fVal / (value_type)1e3; }
value_type Rnd(value_type v) { return v*std::rand()/(value_type)(RAND_MAX+1.0); }
value_type Not(value_type v) { return v==0; }
value_type Add(value_type v1, value_type v2) { return v1+v2; }
value_type Mul(value_type v1, value_type v2) { return v1*v2; }
value_type StrFun2(const char_type *v1, value_type v2,value_type v3) 
{ 
  mu::console() << v1 << std::endl;
  return v2+v3; 
}

value_type Ping() { mu::console() << "ping\n"; return 0; }

mu::value_type SampleQuery(const char_type *szMsg) 
{
  if (szMsg) 
    mu::console() << szMsg << std::endl;

  return 999;
};

//---------------------------------------------------------------------------
// Factory function for creating new parser variables
// This could as well be a function performing database queries.
double* AddVariable(const char_type *a_szName, void *a_pUserData)
{
  // I don't want dynamic allocation here, so i used this static buffer
  // If you want dynamic allocation you must allocate all variables dynamically
  // in order to delete them later on. Or you find other ways to keep track of 
  // variables that have been created implicitely.
  static double afValBuf[100];  
  static int iVal = 0;          

  mu::console() << _T("Generating new variable \"") 
                << a_szName << _T("\" (slots left: ")
                << 99-iVal << _T(")")
                << _T(" User data pointer is:") 
                << std::hex << a_pUserData <<endl;

  afValBuf[iVal] = 0;
  if (iVal>=99)
    throw mu::ParserError( _T("Variable buffer overflow.") );

  return &afValBuf[iVal++];
}

//---------------------------------------------------------------------------
void Splash()
{
  mu::console() << _T( "                 __________                                       \n");
  mu::console() << _T( "    _____   __ __\\______   \\_____  _______  ______  ____ _______\n");
  mu::console() << _T( "   /     \\ |  |  \\|     ___/\\__  \\ \\_  __ \\/  ___/_/ __ \\\\_  __ \\ \n");
  mu::console() << _T( "  |  Y Y  \\|  |  /|    |     / __ \\_|  | \\/\\___ \\ \\  ___/ |  | \\/ \n");
  mu::console() << _T( "  |__|_|  /|____/ |____|    (____  /|__|  /____  > \\___  >|__|    \n");
  mu::console() << _T( "        \\/                       \\/            \\/      \\/         \n");
  mu::console() << _T( "  (C) 2004-2007 Ingo Berg                                         \n");
  mu::console() << _T( "                                                                  \n");
}

//---------------------------------------------------------------------------
void SelfTest()
{
  mu::console() << _T( "-----------------------------------------------------------\n");
  mu::console() << _T( "Running test suite:\n");
  mu::console() << _T( "-----------------------------------------------------------\n");

  mu::Test::ParserTester pt;
  pt.Run();

  mu::console() << _T( "-----------------------------------------------------------\n");
  mu::console() << _T( "Commands:\n");
  mu::console() << _T( "  list var     - list parser variables\n");
  mu::console() << _T( "  list exprvar - list expression variables\n");
  mu::console() << _T( "  list const   - list all numeric parser constants\n");
  mu::console() << _T( "  exit         - exits the parser\n");
  mu::console() << _T( "Constants:\n");
  mu::console() << _T( "  \"_e\"   2.718281828459045235360287\n");
  mu::console() << _T( "  \"_pi\"  3.141592653589793238462643\n");
  mu::console() << _T( "-----------------------------------------------------------\n");
  mu::console() << _T( "Enter a formula or a command:\n");
}

//---------------------------------------------------------------------------
void ListVar(const mu::ParserBase &parser)
{
  // Query the used variables (must be done after calc)
  mu::varmap_type variables = parser.GetVar();
  if (!variables.size())
    return;

  cout << "\nParser variables:\n";
  cout <<   "-----------------\n";
  cout << "Number: " << (int)variables.size() << "\n";
  varmap_type::const_iterator item = variables.begin();
  for (; item!=variables.end(); ++item)
    mu::console() << _T("Name: ") << item->first << _T("   Address: [0x") << item->second << _T("]\n");
}

//---------------------------------------------------------------------------
void ListConst(const mu::ParserBase &parser)
{
  mu::console() << _T("\nParser constants:\n");
  mu::console() << _T("-----------------\n");

  mu::valmap_type cmap = parser.GetConst();
  if (!cmap.size())
  {
    mu::console() << _T("Expression does not contain constants\n");
  }
  else
  {
    valmap_type::const_iterator item = cmap.begin();
    for (; item!=cmap.end(); ++item)
      mu::console() << _T("  ") << item->first << _T(" =  ") << item->second << _T("\n");
  }
}

//---------------------------------------------------------------------------
void ListExprVar(const mu::ParserBase &parser)
{
  string_type sExpr = parser.GetExpr();
  if (sExpr.length()==0)
  {
    cout << _T("Expression string is empty\n");
    return;
  }

  // Query the used variables (must be done after calc)
  mu::console() << _T("\nExpression variables:\n");
  mu::console() <<   _T("---------------------\n");
  mu::console() << _T("Expression: ") << parser.GetExpr() << _T("\n");

  varmap_type variables = parser.GetUsedVar();
  if (!variables.size())
  {
    mu::console() << _T("Expression does not contain variables\n");
  }
  else
  {
    mu::console() << _T("Number: ") << (int)variables.size() << _T("\n");
    mu::varmap_type::const_iterator item = variables.begin();
    for (; item!=variables.end(); ++item)
      mu::console() << _T("Name: ") << item->first << _T("   Address: [0x") << item->second << _T("]\n");
  }
}

//---------------------------------------------------------------------------
/** \brief Check for external keywords.
*/
bool CheckKeywords(const mu::char_type *a_szLine, mu::ParserBase &a_Parser)
{
  string_type sLine(a_szLine);

  if ( sLine == _T("quit") )
  {
    exit(0);
  }
  else if ( sLine == _T("list var") )
  {
    ListVar(a_Parser);
    return true;
  }
  else if ( sLine == _T("list const") )
  {
    ListConst(a_Parser);
    return true;
  }
  else if ( sLine == _T("list exprvar") )
  {
    ListExprVar(a_Parser);
    return true;
  }
  else if ( sLine==_T("list const") )
  {
    ListConst(a_Parser);
    return true;
  }

  return false;
}

//---------------------------------------------------------------------------
void Calc()
{
  mu::Parser  parser;
  mu::ParserInt int_parser;

  // Add some variables
  value_type  vVarVal[] = { 1, 2 }; // Values of the parser variables
  parser.DefineVar( _T("a"), &vVarVal[0]);  // Assign Variable names and bind them to the C++ variables
  parser.DefineVar( _T("b"), &vVarVal[1]);
  parser.DefineStrConst( _T("strBuf"), _T("hello world") );

  // Add user defined unary operators
  parser.DefinePostfixOprt( _T("M"), Mega);
  parser.DefinePostfixOprt( _T("m"), Milli);
  parser.DefineInfixOprt( _T("!"), Not);
  parser.DefineFun( _T("query"), SampleQuery, false);
  parser.DefineFun( _T("rnd"), Rnd, false);     // Add an unoptimizeable function
  parser.DefineFun( _T("strfun2"), StrFun2, false); // Add an unoptimizeable function
  parser.DefineFun( _T("ping"), Ping, false);

  parser.DefineOprt( _T("add"), Add, 0);
  parser.DefineOprt( _T("mul"), Mul, 1);
  parser.DefineOprt( _T("$"), Mul, 1);

  // Define the variable factory
  parser.SetVarFactory(AddVariable, &parser);

  for(;;)
  {
    try
    {
      string_type sLine;
      std::getline(mu::console_in(), sLine);

      if (CheckKeywords(sLine.c_str(), parser)) 
        continue;

//#define MUP_EXAMPLE_INT_PARSER
#ifdef MUP_EXAMPLE_INT_PARSER
      int_parser.SetExpr(sLine);
      mu::console() << int_parser.Eval() << "\n";
#else
      parser.SetExpr(sLine);
      parser.Eval();
      mu::console() << parser.Eval() << "\n";
#endif
    }
    catch(mu::Parser::exception_type &e)
    {
      mu::console() << _T("\nError:\n");
      mu::console() << _T("------\n");
      mu::console() << _T("Message:     ")   << e.GetMsg()   << _T("\n");
      mu::console() << _T("Expression:  \"") << e.GetExpr()  << _T("\"\n");
      mu::console() << _T("Token:       \"") << e.GetToken()    << _T("\"\n");
      mu::console() << _T("Position:    ")   << (int)e.GetPos() << _T("\n");
      mu::console() << _T("Errc:        ")   << std::dec << e.GetCode() << _T("\n");
    }
  } // while running
}

//---------------------------------------------------------------------------
int main(int, char**)
{
  using namespace mu;

  Splash();
  SelfTest();

  try
  {
    Calc();
  }
  catch(Parser::exception_type &e)
  {
    // Only erros raised during the initialization will end up here
    // formula related errors are treated in Calc()
    console() << _T("Initialization error:  ") << e.GetMsg() << endl;

    string_type sBuf;
    console_in() >> sBuf;
  }

  return 0;
}
