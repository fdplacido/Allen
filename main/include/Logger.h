#pragma once

#define verbose_cout logger::logger(logger::verbose)
#define debug_cout logger::logger(logger::debug)
#define info_cout logger::logger(logger::info)
#define warning_cout logger::logger(logger::warning)
#define error_cout logger::logger(logger::error)

#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include "LoggerCommon.h"

// Dumb type, just making constructor public
class FileStdLogger : public std::ostream {
public:
  FileStdLogger() : std::ostream() {}
};

// This is relatively simple
class MessageLogger : public std::streambuf {
public:
  std::string _buf;
  std::ofstream* _file_io = nullptr;
  FileStdLogger* _file_std_io = nullptr;
  std::streambuf* _old = nullptr;

  MessageLogger(std::ofstream* file_io, FileStdLogger* file_std_io) :
    _buf(""), _file_io(file_io), _file_std_io(file_std_io)
  {
    // Override the previous read buffer
    _old = _file_std_io->rdbuf(this);
  }

  ~MessageLogger() { _file_std_io->rdbuf(_old); }

  int overflow(int c) override
  {
    if (c == '\n') {
      std::cout << _buf << std::endl;
      (*_file_io) << _buf << std::endl;
      _buf = "";
    }
    else
      _buf += c;

    return c;
  }
};

class VoidLogger : public std::streambuf {
public:
  std::ostream* _void_stream = nullptr;
  std::streambuf* _old = nullptr;

  VoidLogger(std::ostream* void_stream) : _void_stream(void_stream) { _old = _void_stream->rdbuf(this); }

  ~VoidLogger() { _void_stream->rdbuf(_old); }

  int overflow(int c) override
  {
    // Just don't do anything
    return c;
  }
};

namespace logger {
  class Logger {
  public:
    int verbosityLevel = 3;
    FileStdLogger discardStream;
    VoidLogger* discardLogger = nullptr;
    Logger() { discardLogger = new VoidLogger(&discardStream); }
  };

  std::ostream& logger(int requestedLogLevel);

  extern Logger ll;
} // namespace logger
