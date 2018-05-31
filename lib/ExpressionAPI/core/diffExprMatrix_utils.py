import os
import uuid
import re
import json
import numpy
from pprint import pprint, pformat

from Workspace.WorkspaceClient import Workspace
from GenomeAnnotationAPI.GenomeAnnotationAPIClient import GenomeAnnotationAPI

class DiffExprMatrixUtils:
    """
     Constains a set of functions for expression levels calculations.
    """

    PARAM_IN_WS_NAME = 'workspace_name'
    PARAM_IN_OBJ_NAME = 'output_obj_name'
    PARAM_IN_DIFFEXPMATSET_REF = 'diffExprMatrixSet_ref'

    def __init__(self, config, logger=None):
        self.config = config
        self.logger = logger
        self.scratch = os.path.join(config['scratch'], 'DEM_' + str(uuid.uuid4()))
        self.ws_url = config['workspace-url']
        self.ws_client = Workspace(self.ws_url)
        self.serviceWizardURL = config['srv-wiz-url']
        self._mkdir_p(self.scratch)
        pass

    def _mkdir_p(self, path):
        """
        _mkdir_p: make directory for given path
        """
        if not path:
            return
        try:
            os.makedirs(path)
        except OSError as exc:
            if exc.errno == errno.EEXIST and os.path.isdir(path):
                pass
            else:
                raise

    def process_params(self, params):
        """
        validates params passed to gen expression matrix method
        """
        for p in [self.PARAM_IN_DIFFEXPMATSET_REF]:
            if p not in params:
                raise ValueError('"{}" parameter is required, but missing'.format(p))

    def get_expressionset_data(self, expressionset_ref):

        expr_set_obj = self.ws_client.get_objects2(
            {'objects': [{'ref': expressionset_ref}]})['data'][0]

        expr_set_obj_type = expr_set_obj.get('info')[2]
        expr_set_data = dict()
        expr_set_data['ws_name'] = expr_set_obj.get('info')[7]
        expr_set_data['obj_name'] = expr_set_obj.get('info')[1]

        if re.match('KBaseRNASeq.RNASeqExpressionSet-\d.\d', expr_set_obj_type):
            expr_set_data['genome_ref'] = expr_set_obj['data']['genome_id']
            expr_obj_refs = list()
            for expr_obj in expr_set_obj['data']['mapped_expression_ids']:
                expr_obj_refs.append(expr_obj.values()[0])
            expr_set_data['expr_obj_refs'] = expr_obj_refs

        elif re.match('KBaseSets.ExpressionSet-\d.\d', expr_set_obj_type):
            items = expr_set_obj.get('data').get('items')
            expr_obj_refs = list()
            for item in items:
                expr_obj_refs.append(item['ref'])
            expr_obj = self.ws_client.get_objects2(
                {'objects': [{'ref': expr_obj_refs[0]}]})['data'][0]
            expr_set_data['genome_ref'] = expr_obj['data']['genome_id']
            expr_set_data['expr_obj_refs'] = expr_obj_refs
        else:
            raise TypeError(self.PARAM_IN_EXPSET_REF + ' should be of type ' +
                            'KBaseRNASeq.RNASeqExpressionSet ' +
                            'or KBaseSets.ExpressionSet')
        return expr_set_data

    def get_diffexpr_matrixset(self, params, token):

        self.ws_client = Workspace(self.ws_url, token=token)

        col_names = {'gene_id': 'gene',
                     'log2_fold_change': 'log2fc_f',
                     'p_value': 'p_value_f',
                     'q_value': 'q_value'}

        json_fields = ['log2fc_f', 'p_value_f', 'q_value']

        self.process_params(params)

        diffexprmatset_list = list()
        diffexprmatset_ref = params.get(self.PARAM_IN_DIFFEXPMATSET_REF)

        diffexprmatset_obj = self.ws_client.get_objects2(
                                {'objects': [{'ref': diffexprmatset_ref}]})['data'][0]

        items = diffexprmatset_obj.get('data').get('items')
        diffexprmat_refs = list()

        for item in items:
            diffexprmat_refs.append(item['ref'])
            self.logger.info('DiffExprMatrix ref: ' + item['ref'])

        for diffexprmat_ref in diffexprmat_refs:
            diffexprmat_dict = dict()
            diffexprmat_obj = self.ws_client.get_objects2(
                                {'objects': [{'ref': diffexprmat_ref}]})['data'][0]
            diffexprmat = diffexprmat_obj.get('data')
            diffexprmat_dict['condition_1'] = diffexprmat.get('condition_mapping').keys()[0]
            diffexprmat_dict['condition_2'] = diffexprmat.get('condition_mapping').values()[0]
            voldata = list()
            data = diffexprmat.get('data')

            for row_index, row_id in enumerate(data.get('row_ids')):
                row_data = dict()
                row_data['gene'] = row_id
                values = data.get('values')[row_index]
                for col_index in range(len(values)):
                    row_data[json_fields[col_index]] = values[col_index]

                voldata.append(row_data)

            diffexprmat_dict['voldata'] = voldata
            diffexprmatset_list.append(diffexprmat_dict)

        return diffexprmatset_list


    def get_matrix_stats( self, raw_row ):
        """
        returns a list of [ min, max, mean, std.dev, is_data_missing] for one row of conditional 
        expression values
        """
        has_missing = "No"
        row = []
        for r in raw_row:
            if r == None or numpy.isnan( r ):     # careful here - r can be 0 which is a legitimate value
                has_missing = "Yes"
            else:
                row.append(r)

        if len( row ) < 1:
            return( [ 'NA', 'NA', 'NA', 'NA', 'Yes' ] )

        if len( row ) == 1:
           sd = 0
        else:
           sd = numpy.std( row, ddof=1 )
        return( [ min( row ), max( row ), numpy.mean( row ), sd, has_missing ] )


    def convert_dem_to_dict( self, dem ):
        """
        returns a dict that maps feature_id -> [ fc, q ]
        """
        row_ids = dem.get( 'row_ids' )
        vals = dem.get( 'values' )
  
        n_rows = len( row_ids )
        if ( len( vals ) != n_rows ):
            raise Exception( "length discrepancy in differential expression matrix: {0} row_ids but {1} values".format( n_rows, len( fvals ) ) )

        dem_dict = {}
        for _id, val in zip(row_ids, vals):
            dem_dict[_id] = [ val[0], val[2] ]  # [fc,q]. (not bothering to check for dups here)

        return dem_dict


    def get_enhancedFEM( self, params, tok ):
        """
        implements get_enhancedFilteredExpressionMatrix() method
        """

        if 'fem_object_ref' not in params:
            raise ValueError( "fem_object_ref parameter not given to get_enhancedFilteredExpressionMatrix" )

        fem_object_ref = params.get( 'fem_object_ref' )

        fem_obj_ret = self.ws_client.get_objects2(
                       {'objects': [{'ref': fem_object_ref }]})['data'][0]
        fem = fem_obj_ret.get( 'data' )
        prov = fem_obj_ret.get( 'provenance')[0]

        # create the enhanced FEM, starting with the FEM

        efem = {}
        for k in [ 'genome_ref', 'scale', 'type' ]:
            efem[k] = fem.get( k )

        efem['data'] = {}
        efem['data']['col_ids'] = [ "description", 
                                    "fold-change",
                                    "q-value",
                                    "min",
                                    "max",
                                    "mean",
                                    "std_dev",
                                    "is_missing_values" ]
        efem['data']['column_labels'] =[ "Description", 
                                         "Fold change",
                                         "Q value",
                                         "Min. expression",
                                         "Max. expression",
                                         "Mean expression",
                                         "Std. dev.",
                                         "Missing values?" ]
        fm = fem.get('data')
        efem['data']['row_ids'] = fm.get('row_ids')
        efem['data']['values' ] = []
        n_efem_rows = len( efem['data']['row_ids'] )
        fvals = fm.get('values')
        if ( len( fvals ) != n_efem_rows ):
            raise Exception( "length discrepancy in filtered expression matrix: {0} row_ids but {1} values".format( n_efem_rows, len( fvals ) ) )

        # Get genome object and feature descriptions as a handy feature-indexed dict

        # moved from constructor
        gaa = GenomeAnnotationAPI( self.serviceWizardURL, token=tok )
        feat_dict = gaa.get_feature_functions( { 'ref': fem.get( 'genome_ref' ), 'feature_id_list': None } )

        # if this FEM has a "resolved_ws_objects" record in its provenance,
        # then that should be a list of one DEM reference from which we get the FC and q values
        # as a feature (=row_id) -indexed dict.

        if fem.get( 'diff_expr_matrix_ref' ):
            dem_ref = fem.get( 'diff_expr_matrix_ref' )
            dem_obj_ret = self.ws_client.get_objects2(
                       {'objects': [{'ref': dem_ref }]})['data'][0]
      
            dem = dem_obj_ret.get( 'data' )
            dem_dict = self.convert_dem_to_dict( dem.get('data') )  # convert to dictionary for quick lookups
        else:
            dem_dict = {}   # empty dictionary

        # for each row

        for row_id, fm_val_row in zip( fm.get('row_ids'), fvals ):

            # make a new row with NA for description, FC and q

            new_values_row =  [ 'NA', 'NA', 'NA' ] + self.get_matrix_stats( fm_val_row )

            # if we have a description for this feature (row_id) put it in the first column

            desc = feat_dict.get( row_id )
            if desc:
                new_values_row[0] = desc     # leave as 'NA' if no entry in feat_dict

            # if we have a DEM entry for this row, put FC and q into 2nd and 3rd columns
            d = dem_dict.get( row_id )
            if d:
                new_values_row[1], new_values_row[2] = d

            # finally, add this row to the eFEM

            efem['data']['values'].append( new_values_row )


        return efem





