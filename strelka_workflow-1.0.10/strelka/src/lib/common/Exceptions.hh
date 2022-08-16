// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2009-2013 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/downloads/sequencing/licenses/>.
//

/**
 ** \brief Declaration of the common exception mechanism.
 **
 ** All exceptions must carry the same data (independently of the
 ** exception type) to homogenixe the reporting and processing of
 ** errors.
 **
 ** \author Come Raczy
 **/

#pragma once

#include <string>
#include <stdexcept>
#include <ios>
#include <boost/cerrno.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/exception/all.hpp>
#include <boost/throw_exception.hpp>

namespace illumina
{
namespace common
{

/**
 ** \brief Virtual base class to all the exception classes in CASAVA.
 **
 ** Use BOOST_THROW_EXCEPTION to get the contect info (file, function, line)
 ** at the throw site.
 **/
class ExceptionData : public boost::exception
{
public:
    ExceptionData(int errorNumber=0, const std::string& message="");
    ExceptionData(const ExceptionData& e) : boost::exception(e), errorNumber_(e.errorNumber_), message_(e.message_) {}
    virtual ~ExceptionData() throw ()
    {
    }
    int getErrorNumber() const
    {
        return errorNumber_;
    }
    std::string getMessage() const
    {
        return message_;
    }
    std::string getContext() const;
private:
    const int errorNumber_;
    const std::string message_;
    ExceptionData& operator=(const ExceptionData&);
};

class IlluminaException: public std::exception, public ExceptionData
{
public:
    IlluminaException(int errorNumber, const std::string& message) : ExceptionData(errorNumber, message) {}
    IlluminaException(const IlluminaException& e) : std::exception(e), ExceptionData(e) {}
private:
    IlluminaException& operator=(const IlluminaException&);
};

class IoException: public std::ios_base::failure, public ExceptionData
{
public:
    IoException(int errorNumber, const std::string& message);
};

/**
 ** \brief Exception thrown when the client supplied and unsupported version number.
 **
 ** Particularly relevant to data format and software versions
 ** (Pipeline, IPAR, Phoenix, etc.). It should not be used in
 ** situations where the client didn't have the possibility to check
 ** the version (for instance when reading the version of a data
 ** format from the header of a file).
 **
 **/
class UnsupportedVersionException: public std::logic_error, public ExceptionData
{
public:
    UnsupportedVersionException(const std::string& message);
};

/**
 ** \brief Exception thrown when the client supplied an invalid parameter.
 **
 **/
class InvalidParameterException: public std::logic_error, public ExceptionData
{
public:
    InvalidParameterException(const std::string& message);
};

/**
 ** \brief Exception thrown when an invalid command line option was detected.
 **
 **/
class InvalidOptionException: public std::logic_error, public ExceptionData
{
public:
    InvalidOptionException(const std::string& message);
};

/**
 ** \brief Exception thrown when a method invocation violates the pre-conditions.
 **
 **/
class PreConditionException: public std::logic_error, public ExceptionData
{
public:
    PreConditionException(const std::string& message);
};

/**
 ** \brief Exception thrown when a method invocation violates the post-conditions.
 **
 **/
class PostConditionException: public std::logic_error, public ExceptionData
{
public:
    PostConditionException(const std::string& message);
};

}
}
