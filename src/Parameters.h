/* -------------------------------------------------------------------------- *
 *                                  QCDyn                                     *
 * -------------------------------------------------------------------------- *
 * This is part of the Quantum Classical Dynamics (QCDyn) program.            *
 *                                                                            *
 * Author: Zhubin Hu @Sun Group @NYU-SH                                       *
 * Last updated: Nov. 26, 2021                                                *
 * -------------------------------------------------------------------------- */

#pragma once
#include "Tools.h"

/**
 * This singleton class is used to store global parameters.
 * You can use Parameters::get() to construct a unique object.
 * The function: loadParamters() can be used to load parameters from input control
 * file or command line.
 * The functions: getStr(), getInt(), getDouble() and getBool() can be used to get
 * the value with specified type that you want.
 */
class Parameters {
public:
    /**
     * Construct a static Parameters object.
     */
    static Parameters& get() {
        static Parameters instance;
        return instance;
    }
    /**
     * Load parameters from input control file and command line.
     * This function will call loadParamFromFile() and loadParamFromString() to
     * load paramters.
     *
     * If an input control file in which contains parameters is provided in command
     * line (such as control.inp), load parameters from the input file firstly which
     * will add new parameters (isNew=true) or overwrite the default values (isNew=true)
     * depends on the tag called isNew passed to loadParamFromFile().
     * Then load the parameters provided in the command line (such as --job=NVE),
     * which will overwrite the existing value.
     *
     * If no argument is provided, it will print the usage and default paramters
     * (if not empty), then exit.
     *
     * @param argc   the number of arguments from command line
     * @param argv   the arguments from command line, argc[0] is program itself
     * @param isNew  true: add new parameters; false: set values for existing parmeters
     */
    void loadParamters(int argc, char *argv[], bool isNew = false);
    /**
     * Add new parameters or set values for existing parmeters from control file.
     *
     * File format:
     * 1st column: key; 2nd column: value; 3rd column (optional): type;
     * Each column should be divided with at least a space.
     * No space used in the key, value, and type.
     * A line starting with "#" or empty will be ignored.
     *
     * If isNew = true, then the type of parameter must be provided, which means
     * initialize the Parameters class members from control file.
     * If isNew = false (default), the type of parameter is not necessary, which
     * should be used when the Parameters class members have been initialized in
     * within the program.
     *
     * @param file   file name of input control file
     * @param isNew  true: add new parameters; false: set values for existing parmeters
     */
    void loadParamFromFile(const std::string& file, bool isNew = false);
    /**
     * Set values for existing parmeters from string with a style like "--id=123",
     * which can be used to load parameters from command line.
     *
     * If no such key exists, throw an exception.
     *
     * @param command   the string like "--id=123" including key and value
     */
    void loadParamFromString(const std::string& command);
    /**
     * Add a new paramter entry.
     *
     * The parameter of type is string and will be converted according to the given type.
     * The allowed types are "str", "int", "double", "bool".
     * For a bool type, the allowed values are "true" and "false", or 1 and 0.
     *
     * If the key already exists or unsupported type is provided, throw an exception.
     *
     * @param key    the key of parameter
     * @param value  the value for the key
     * @param type   the type for the key
     */
    void addParam(const std::string& key, const std::string& value, const std::string& type);
    /**
     * Set value for a existing parmeter entry with a set of key and value.
     *
     * If no such key exists, throw an exception.
     *
     * @param key    the key of parameter to be set
     * @param value  the value for the key
     */
    void setValue(const std::string& key, const std::string& value);
    /**
     * Print the information of one parameter with a given key.
     *
     * If no such key exists, throw an exception.
     *
     * @param key    the key of parameter to be printed
     */
    void printParam(const std::string& key) const;
    /**
     * Print all the parameters with current value in the adding order.
     */
    void printAllParam() const;
    /**
     * Get the type of a given key.
     *
     * We used the .at(key) member function of map here, which will check the existence
     * of key, i.e., if no such element exists, throw an exception. Note that, when
     * using the operator "[key]", if no such element exists, it will insert specified
     * element which will modify the map but without value.
     *
     * If no such key exists, throw an exception.
     *
     * @param key    the key of parameter
     * @return       the type of this parameter
     */
    std::string getType(const std::string& key) const;
    /**
     * Get the string value of a given key.
     * If no such key exists, throw an exception.
     *
     * @param key    the key of parameter to get value
     * @return       the string value of this parameter
     */
    std::string getStr(const std::string& key) const;
    /**
     * Get the integer value of a given key.
     * If no such key exists, throw an exception.
     *
     * @param key    the key of parameter to get value
     * @return       the integer value of this parameter
     */
    int getInt(const std::string& key) const;
    /**
     * Get the double value of a given key.
     * If no such key exists, throw an exception.
     *
     * @param key    the key of parameter to get value
     * @return       the double value of this parameter
     */
    double getDouble(const std::string& key) const;
    /**
     * Get the bool value of a given key.
     * If no such key exists, throw an exception.
     *
     * @param key    the key of parameter to get value
     * @return       the bool value of this parameter
     */
    bool getBool(const std::string& key) const;

private:
    /**
     * Internal class used to create a single entry of parameter.
     *
     * For each entry, there are four values with different types.
     * For a string type, int/double_value is 0, and if string value is "true", the
     * bool_value is true, otherwise, it is false.
     * For a int/double type, if the value is 0, the bool_value is false, otherwise,
     * it is true.
     * For a bool type, when the value is "true", or non-zero number, the bool_value
     * is true, otherwise it is false. The string value will be "true/false", and the
     * int/double value will be 1/0.
     * @private
     */
    class ParamEntry {
    public:
        std::string name;
        std::string type;
        std::string str_value;
        int         int_value;
        double      double_value;
        bool        bool_value;
        /**
         * Construct a paramter entry.
         *
         * The parameter of type is string and will be converted according to the given type.
         * The allowed types are "str", "int", "double", "bool".
         * For a bool type, the allowed values are "true" and "false", or 1 and 0.
         *
         * @param key    the key of parameter
         * @param value  the value for the key
         * @param type   the type for the key
         */
        ParamEntry(const std::string& key, const std::string& value, const std::string& type);
        ParamEntry() {}
        ~ParamEntry() {}
    };
    /**
     * The map is used to store the keys and values of parameters.
     * The vector is used to stores the keys only with the adding order of parameters.
     */
    std::map<std::string, ParamEntry> parameters;
    std::vector<std::string> keys;
    Parameters() {}
    ~Parameters() {}
};

/******************************************************************************
 * The following functions are designed for this code only.                   *
 ******************************************************************************/
/**
 * Initialize the Parameter object used for simulation.
 *
 * @param argc     the number of arguments from command line
 * @param argv     the arguments from command line, argc[0] is program itself
 * @param param    the object used to store global parameters
 */
void InitializeParameters(int argc, char *argv[], Parameters& param);
/**
 * Add all supported keys with the default values and types to Parameters.
 *
 * @param param    the object used to store global parameters
 */
static void AddDefaultParameters(Parameters& param);
/**
 * Set default values of some parameters according to job type.
 *
 * @param param    the object used to store global parameters
 */
static void SetJobParameters(Parameters& param);