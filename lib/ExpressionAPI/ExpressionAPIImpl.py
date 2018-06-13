# -*- coding: utf-8 -*-
#BEGIN_HEADER
import logging
import os
import sys
import time
from core.diffExprMatrix_utils import DiffExprMatrixUtils
#END_HEADER


class ExpressionAPI:
    '''
    Module Name:
    ExpressionAPI

    Module Description:
    
    '''

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "0.1.0"
    GIT_URL = "git@github.com:sean-mccorkle/ExpressionAPI.git"
    GIT_COMMIT_HASH = "c59070bcb64009d699e6c8b14f374152617a299e"

    #BEGIN_CLASS_HEADER
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.__LOGGER = logging.getLogger('ExpressionUtils')
        self.__LOGGER.setLevel(logging.INFO)
        streamHandler = logging.StreamHandler(sys.stdout)
        formatter = logging.Formatter(
            "%(asctime)s - %(filename)s - %(lineno)d - %(levelname)s - %(message)s")
        formatter.converter = time.gmtime
        streamHandler.setFormatter(formatter)
        self.__LOGGER.addHandler(streamHandler)
        self.__LOGGER.info("Logger was set")
        self.config = config
        self.scratch = config['scratch']
        self.ws_url = config['workspace-url']
        self.diffexpr_matrix_utils = DiffExprMatrixUtils(config, self.__LOGGER)
        #END_CONSTRUCTOR
        pass


    def get_differentialExpressionMatrixSet(self, ctx, params):
        """
        :param params: instance of type "getDiffExprMatrixParams" (*
           Following are the required input parameters to get Differential
           Expression Matrix json object *) -> structure: parameter
           "diffExprMatrixSet_ref" of String
        :returns: instance of type "getDiffExprMatrixOutput" -> structure:
           parameter "volcano_plot_data" of unspecified object, parameter
           "json_filepath" of String
        """
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN get_differentialExpressionMatrixSet

        plot_data = self.diffexpr_matrix_utils.get_diffexpr_matrixset(params, ctx['token'])

        returnVal = {'volcano_plot_data': plot_data}

        #END get_differentialExpressionMatrixSet

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method get_differentialExpressionMatrixSet return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]

    def get_enhancedFilteredExpressionMatrix(self, ctx, params):
        """
        :param params: instance of type "getEnhancedFEMParams" (* Input
           parameters and method for getting the enhanced Filtered Expresion
           Matrix for viewing *) -> structure: parameter "fem_object_ref" of
           String
        :returns: instance of type "getEnhancedFEMOutput" -> structure:
           parameter "enhanced_FEM" of type "ExpressionMatrix" (A wrapper
           around a FloatMatrix2D designed for simple matricies of Expression
           data.  Rows map to features, and columns map to conditions.  The
           data type includes some information about normalization factors
           and contains mappings from row ids to features and col ids to
           conditions. description - short optional description of the
           dataset type - ? level, ratio, log-ratio scale - ? probably: raw,
           ln, log2, log10 col_normalization - mean_center, median_center,
           mode_center, zscore row_normalization - mean_center,
           median_center, mode_center, zscore feature_mapping - map from
           row_id to feature id in the genome data - contains values for
           (feature,condition) pairs, where features correspond to rows and
           conditions are columns (ie data.values[feature][condition])
           diff_expr_matrix_ref - added to connect filtered expression matrix
           to differential expression matrix used for filtering @optional
           description row_normalization col_normalization @optional
           genome_ref feature_mapping conditionset_ref condition_mapping
           report diff_expr_matrix_ref @metadata ws type @metadata ws scale
           @metadata ws row_normalization @metadata ws col_normalization
           @metadata ws genome_ref as Genome @metadata ws conditionset_ref as
           ConditionSet @metadata ws length(data.row_ids) as feature_count
           @metadata ws length(data.col_ids) as condition_count) ->
           structure: parameter "description" of String, parameter "type" of
           String, parameter "scale" of String, parameter "row_normalization"
           of String, parameter "col_normalization" of String, parameter
           "genome_ref" of type "ws_genome_id" (The workspace ID for a Genome
           data object. @id ws KBaseGenomes.Genome), parameter
           "feature_mapping" of mapping from String to String, parameter
           "conditionset_ref" of type "ws_conditionset_id" (The workspace ID
           for a ConditionSet data object (Note: ConditionSet objects do not
           yet exist - this is for now used as a placeholder). @id ws
           KBaseExperiments.ConditionSet), parameter "condition_mapping" of
           mapping from String to String, parameter "diff_expr_matrix_ref" of
           String, parameter "data" of type "FloatMatrix2D" (A simple 2D
           matrix of floating point numbers with labels/ids for rows and
           columns.  The matrix is stored as a list of lists, with the outer
           list containing rows, and the inner lists containing values for
           each column of that row.  Row/Col ids should be unique. row_ids -
           unique ids for rows. col_ids - unique ids for columns. values -
           two dimensional array indexed as: values[row][col] @metadata ws
           length(row_ids) as n_rows @metadata ws length(col_ids) as n_cols)
           -> structure: parameter "row_ids" of list of String, parameter
           "col_ids" of list of String, parameter "values" of list of list of
           Double, parameter "report" of type "AnalysisReport" (A basic
           report object used for a variety of cases to mark informational
           messages, warnings, and errors related to processing or quality
           control checks of raw data.) -> structure: parameter
           "checkTypeDetected" of String, parameter "checkUsed" of String,
           parameter "checkDescriptions" of list of String, parameter
           "checkResults" of list of type "boolean" (Indicates true or false
           values, false = 0, true = 1 @range [0,1]), parameter "messages" of
           list of String, parameter "warnings" of list of String, parameter
           "errors" of list of String
        """
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN get_enhancedFilteredExpressionMatrix

        efem = self.diffexpr_matrix_utils.get_enhancedFEM( params, ctx['token'] )

        returnVal = { 'enhanced_FEM': efem }

        #END get_enhancedFilteredExpressionMatrix

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method get_enhancedFilteredExpressionMatrix return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]
    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK",
                     'message': "",
                     'version': self.VERSION,
                     'git_url': self.GIT_URL,
                     'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
