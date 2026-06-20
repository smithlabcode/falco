/* MIT License
 *
 * Copyright (c) 2026 Andrew D Smith
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include "html.hpp"

#include "falco_grade.hpp"
#include "falco_utils.hpp"

#include <config.h>

#define FMT_HEADER_ONLY
#include "fmt/chrono.h"
#include "fmt/format.h"
#include "fmt/ranges.h"

#include <chrono>
#include <string>
#include <vector>

static constexpr auto style = R"(<style type="text/css">
  @media screen {
    div.summary {
      width: 18em;
      position:fixed;
      top: 4em;
      margin:1em 0 0 1em;
  }
  div.main {
    display:block;
    position:absolute;
    overflow:auto;
    height:auto;
    width:auto;
    top:4.5em;
    bottom:2.3em;
    left:18em;
    right:0;
    border-left: 1px solid #CCC;
    padding:0 0 0 1em;
    background-color: white;
    z-index:1;
  }
  div.header {
    background-color: #EEE;
    border:0;
    margin:0;
    padding: 0.2em;
    font-size: 200%;
    position:fixed;
    width:100%;
    top:0;
    left:0;
    z-index:2;
  }
  div.footer {
    background-color: #EEE;
    border:0;
    margin:0;
    padding:0.5em;
    height: 2.5em;
    overflow:hidden;
    font-size: 100%;
    position:fixed;
    bottom:0;
    width:100%;
    z-index:2;
  }
  img.indented {
    margin-left: 3em;
  }
}
@media print {
  img {
    max-width:100% !important;
    page-break-inside: avoid;
  }
  h2, h3 {
    page-break-after: avoid;
  }
  div.header {
    background-color: #FFF;
  }
}
body {
  color: #000;
  background-color: #FFF;
  border: 0;
  margin: 0;
  padding: 0;
}
div.header {
  border:0;
  margin:0;
  padding: 0.5em;
  font-size: 200%;
  width:100%;
}
#header_title {
  display:inline-block;
  float:left;
  clear:left;
}
#header_filename {
  display:inline-block;
  float:right;
  clear:right;
  font-size: 50%;
  margin-right:2em;
  text-align: right;
}
div.header h3 {
  font-size: 50%;
  margin-bottom: 0;
}
div.summary ul {
  padding-left:0;
  list-style-type:none;
}
div.summary ul li img {
  margin-bottom:-0.5em;
  margin-top:0.5em;
}
div.main {
  background-color: white;
}
div.module {
  padding-bottom:3em;
  padding-top:3em;
  border-bottom: 1px solid #990000
}
div.footer {
  background-color: #EEE;
  border:0;
  margin:0;
  padding: 0.5em;
  font-size: 100%;
  width:100%;
}
h2 {
  color: #2a5e8c;
  padding-bottom: 0;
  margin-bottom: 0;
  clear:left;
}
table {
  margin-left: 3em;
  text-align: center;
}
th {
  text-align: center;
  background-color: #000080;
  color: #FFF;
  padding: 0.4em;
}
td {
  font-family: monospace;
  text-align: left;
  background-color: #EEE;
  color: #000;
  padding: 0.4em;
}
img {
  padding-top: 0;
  margin-top: 0;
  border-top: 0;
}
p {
  padding-top: 0;
  margin-top: 0;
}
.pass {
  color : #009900;
}
.warn {
  color : #999900;
}
.fail {
  color : #990000;
}
</style>)";

[[nodiscard]] auto
get_summary(const analysis_grades &ag) -> std::string {
  static constexpr auto summary = R"(
<div class="summary"><h2>Summary</h2>
<ul>
{}
</ul>
</div>
)";
  static constexpr auto list_item =
    R"(<li><a class="{grade}" href="#{module_id}">{module_name}</a></li>)";
  std::vector<std::string> sections;
  for (const auto &name : analysis_grades::names)
    if (ag.is_configured(name))
      sections.emplace_back(fmt::format(
        list_item, fmt::arg("grade", ag.grade(name)),
        fmt::arg("module_id", name), fmt::arg("module_name", ag.label(name))));
  return fmt::format(summary, fmt::join(sections, "\n"));
}

[[nodiscard]] auto
get_html_module(const std::string &name, const std::string &label,
                const std::string &grade,
                const std::string &text) -> std::string {
  static constexpr auto module_template =
    R"(<div class="module">
<h2 class="{module_grade}" id="{module_id}">{module_name}: {module_grade}</h2>
{module_text}
</div>)";
  return fmt::format(module_template,                  //
                     fmt::arg("module_grade", grade),  //
                     fmt::arg("module_id", name),
                     fmt::arg("module_name", label),
                     fmt::arg("module_text", text));
}

static constexpr auto falco_html_body =
  R"(<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1, shrink-to-fit=no">
<title>{filename} - report</title>
<link rel="stylesheet" href="https://stackpath.bootstrapcdn.com/bootstrap/4.3.1/css/bootstrap.min.css" integrity="sha384-ggOyR0iXCbMQv3Xipma34MD+dH/1fQ784/j6cY/iJTQUOhcWr7x9JvoRxT2MZw1T" crossorigin="anonymous">
<link href="https://stackpath.bootstrapcdn.com/font-awesome/4.7.0/css/font-awesome.min.css" rel="stylesheet" integrity="sha384-wvfXpqpZZVQGK6TAh5PVlGOfQNHSoD2xbE+QkPxCAFlNEevoEH3Sl0sibVcOQVnN" crossorigin="anonymous">
{style}
<script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
</head>
<body>
<div class="header">
<div id="header_title">Report</div>
<div id="header_filename">{date}<br/>{filename}</div>
</div>
{summary}
<div class="main">
{modules}
</div>
<div class="footer">Falco {version}</div>
</body>
)";

[[nodiscard]] auto
falco_get_html(const file_info &info, const analysis_grades &grades,
               const std::string &analysis_modules) -> std::string {
  return fmt::format(falco_html_body,
                     fmt::arg("date", std::chrono::system_clock::now()),
                     fmt::arg("filename", info.name),           //
                     fmt::arg("style", style),                  //
                     fmt::arg("summary", get_summary(grades)),  //
                     fmt::arg("modules", analysis_modules),     //
                     fmt::arg("version", VERSION));
}
