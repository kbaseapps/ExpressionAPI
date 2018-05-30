
package us.kbase.kbasefeaturevalues;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import javax.annotation.Generated;
import com.fasterxml.jackson.annotation.JsonAnyGetter;
import com.fasterxml.jackson.annotation.JsonAnySetter;
import com.fasterxml.jackson.annotation.JsonInclude;
import com.fasterxml.jackson.annotation.JsonProperty;
import com.fasterxml.jackson.annotation.JsonPropertyOrder;


/**
 * <p>Original spec-file type: FloatMatrix2D</p>
 * <pre>
 * A simple 2D matrix of floating point numbers with labels/ids for rows and
 * columns.  The matrix is stored as a list of lists, with the outer list
 * containing rows, and the inner lists containing values for each column of
 * that row.  Row/Col ids should be unique.
 * row_ids - unique ids for rows.
 * col_ids - unique ids for columns.
 * values - two dimensional array indexed as: values[row][col]
 * @metadata ws length(row_ids) as n_rows
 * @metadata ws length(col_ids) as n_cols
 * </pre>
 * 
 */
@JsonInclude(JsonInclude.Include.NON_NULL)
@Generated("com.googlecode.jsonschema2pojo")
@JsonPropertyOrder({
    "row_ids",
    "col_ids",
    "values"
})
public class FloatMatrix2D {

    @JsonProperty("row_ids")
    private List<String> rowIds;
    @JsonProperty("col_ids")
    private List<String> colIds;
    @JsonProperty("values")
    private List<List<Double>> values;
    private Map<java.lang.String, Object> additionalProperties = new HashMap<java.lang.String, Object>();

    @JsonProperty("row_ids")
    public List<String> getRowIds() {
        return rowIds;
    }

    @JsonProperty("row_ids")
    public void setRowIds(List<String> rowIds) {
        this.rowIds = rowIds;
    }

    public FloatMatrix2D withRowIds(List<String> rowIds) {
        this.rowIds = rowIds;
        return this;
    }

    @JsonProperty("col_ids")
    public List<String> getColIds() {
        return colIds;
    }

    @JsonProperty("col_ids")
    public void setColIds(List<String> colIds) {
        this.colIds = colIds;
    }

    public FloatMatrix2D withColIds(List<String> colIds) {
        this.colIds = colIds;
        return this;
    }

    @JsonProperty("values")
    public List<List<Double>> getValues() {
        return values;
    }

    @JsonProperty("values")
    public void setValues(List<List<Double>> values) {
        this.values = values;
    }

    public FloatMatrix2D withValues(List<List<Double>> values) {
        this.values = values;
        return this;
    }

    @JsonAnyGetter
    public Map<java.lang.String, Object> getAdditionalProperties() {
        return this.additionalProperties;
    }

    @JsonAnySetter
    public void setAdditionalProperties(java.lang.String name, Object value) {
        this.additionalProperties.put(name, value);
    }

    @Override
    public java.lang.String toString() {
        return ((((((((("FloatMatrix2D"+" [rowIds=")+ rowIds)+", colIds=")+ colIds)+", values=")+ values)+", additionalProperties=")+ additionalProperties)+"]");
    }

}
