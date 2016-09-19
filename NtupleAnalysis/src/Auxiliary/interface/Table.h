#ifndef Table_h
#define Table_h

#include <map>
#include <vector> 
#include <algorithm>
#include <numeric>
#include "constants.h"

using namespace std;
using namespace constants;

// Type definitions
typedef map<int, string> m_Row_To_String;
typedef map<int, string>::iterator it_Row_To_String;
typedef map<int, vector<int> > m_Row_To_Rows;
typedef map<int, m_Row_To_String> m_RowColumn_To_String;
typedef map<int, m_Row_To_String>::iterator it_RowColumn_To_String;


class Table{

 public:
  // Default and overloaded constructors & destructors
  Table() {}; 
  Table(string titleRow, string format, string tableSpecs="");
  ~Table() {};

  // Member Functions
  void SetColumnForRow(int iRow, int iColumn, string newText){ m_RowColumnToString[iRow][iColumn] = newText; }
  string GetRow(int iRow);
  string GetColumnForRow(int iRow, int iColumn){ return m_RowColumnToString[iRow][iColumn]; }
  string GetColumnInRow(string rowText, int iColumn);
  string GetMergedColumnsInRow(int iRow, int iColumnStart, int iColumnEnd);
  string GetTableFormat(void){return format_;}
  string GetTitleRow(void);
  int GetColumnsInString(string title);
  int GetNumberOfColumns(int iRow=-1);
  int GetNumberOfColumnsIncludingDeleted(int iRow=-1);
  int GetNumberOfRows(void);
  int GetNumberOfRowsIncludingDeleted(void);
  bool IsDeletedColumn(int iColumn){ return find(deletedColumns_.begin(), deletedColumns_.end(), iColumn) != deletedColumns_.end(); }
  bool IsDeletedRow(int iRow){ return find(deletedRows_.begin(), deletedRows_.end(), iRow) != deletedRows_.end(); }
  void AddBottomRow(string text){ rowsBottom_.push_back(text); }
  void AddRowColumn(int iRow, string text);
  void AddTopRow(string text){ rowsTop_.push_back(text); }
  void ConvertToFinalState(void);
  void DeleteColumn(int iColumn);
  void DeleteRow(int iRow);
  void InitVars(string titleRow, string format, string tableSpecs="");
  void Print(void);
  void PrintColumn(int iColumn=-1);
  void PrintHorizontalLine(void){ cout << hLine_ << endl; } 
  void PrintRow(int iRow=-1);
  void ReplaceStringInTable(string before, string after);
  void ReplaceStringWithOccurences(int iRowStart, int iRowEnd, int iColumnStart, int iColumnEnd, vector<string> keyWords);
  void SaveToFile(const char *fileName, const char *fileOptions);
  void SetTableWidth(int width){ tableWidth_ = width; } 
    
  
  private:
  // Member Functions
  int _GetColumnWidthForRow(int iRow, int iColumn){ return GetColumnForRow(iRow, iColumn).length(); }
  int _GetNumberOfCharacters(string text);
  string _SetupFormat(string titleRow);
  void _AppendRowEndToEachRow(void);
  void _BeginTabular_(void);
  void _ConvertStringToLatex(string &text);
  void _ConvertTitleRowToLatex(string &titleRow);
  void _DetermineMaxColumnWidths(void);
  void _EndTabular_(void);
  void _FindStringInString(string stringToFind, string stringToSearch);
  void _IsValidColumnIndex(int iColumn);
  void _IsValidRowIndex(int iRow);
  void _IsValidNewRowIndex(int iRow);
  void _PrintAux(void);
  void _PrintRow(vector<string> row){ for (size_t i = 0; i < row.size(); i++){ cout << row.at(i) << endl;} }
  void _PrintTableRows(void);
  void _PrintTitleRow(void);
  void _SetCommentLine(string text){ commentLine_ = text; }
  void _SetDelimiter(string text){ delimiter_ = text; }
  void _SetHorizontalLine(string text){ hLine_ = text; }
  void _SetRowEnd(string text){ rowEnd_ = text; }
  void _SetTableSpecs(string text){ tableSpecs_ = text; }
  void _SetupLatex(void);
  void _SetupText(void);
  
  // Variables  
  bool bCalledPrintAux_;
  int tableWidth_;
  string format_;
  string tableSpecs_;
  string delimiter_;
  string rowEnd_;
  string hLine_;
  string commentLine_;
  vector<string> rowsTop_; 
  vector<string> rowsBottom_;
  vector<int> columnWidths_;
  vector<int> deletedRows_;
  vector<int> deletedColumns_;
  m_Row_To_String m_RowToString;
  m_RowColumn_To_String m_RowColumnToString;

};

#endif
