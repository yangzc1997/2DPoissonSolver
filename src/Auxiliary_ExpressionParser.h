// Auxiliary_ExpressionParser.h
#ifndef AUXILIARY_EXPRESSION_PARSER_H
#define AUXILIARY_EXPRESSION_PARSER_H

#include "UsingAlias.h"
#include <string>

namespace Poisson {

namespace Auxiliary {

// 表达式解析模块
namespace ExpressionParser {
    
    fuxy parse_expression_fuxy(const std::string& expr_str);
    fu parse_expression_fu(const std::string& expr_str);
    fxy parse_expression_fxy(const std::string& expr_str);

} // namespace ExpressionParser

} // namespace Auxiliary

} // namespace Poisson

#endif // AUXILIARY_EXPRESSION_PARSER_H