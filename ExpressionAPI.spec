/*
A KBase module: ExpressionAPI
*/

#include <KBaseFeatureValues.spec>

module ExpressionAPI {
    /*
        Get Differential Expression Matrix from Expression Set input
    */

    /**
        Following are the required input parameters to get Differential Expression Matrix json object
    **/

    typedef structure {

        string      diffExprMatrixSet_ref;

    } getDiffExprMatrixParams;

    typedef structure {

        UnspecifiedObject   volcano_plot_data;
        string              json_filepath;

    } getDiffExprMatrixOutput;

    funcdef  get_differentialExpressionMatrixSet(getDiffExprMatrixParams params)
                                        returns (getDiffExprMatrixOutput)
                                        authentication required;

    /**
        Input parameters and method for getting the enhanced Filtered Expresion Matrix
        for viewing
    **/

    typedef structure {
        string  fem_object_ref;
    } getEnhancedFEMParams;

    typedef structure {
        KBaseFeatureValues.ExpressionMatrix  enhanced_FEM;
    } getEnhancedFEMOutput;


    funcdef  get_enhancedFilteredExpressionMatrix( getEnhancedFEMParams params )
                                   returns (getEnhancedFEMOutput)
                                   authentication required;
                                    

};
