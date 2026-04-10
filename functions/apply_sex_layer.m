function param = apply_sex_layer(param, sex)
%APPLY_SEX_LAYER Select male or female Holmes model file.

    switch lower(sex)
        case 'male'
            param.model = @model_ToRORd_Land_Male;
        case 'female'
            param.model = @model_ToRORd_Land_Female;
        otherwise
            error('Unknown sex: %s', sex);
    end
end