#include "genulens/model/parameters.hpp"

namespace genulens::model {

const ModelParameters &default_model_parameters()
{
    static const ModelParameters params{};
    return params;
}

} // namespace genulens::model

