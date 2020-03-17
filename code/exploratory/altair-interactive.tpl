{% extends 'full.tpl' %}

{% block header %}
  <script src="https://cdn.jsdelivr.net/npm/vega@3"></script>
  <script src="https://cdn.jsdelivr.net/npm/vega-lite@2"></script>
  <script src="https://cdn.jsdelivr.net/npm/vega-embed@3"></script>
  {{super()}}
{% endblock header %}

{%- block data_priority scoped -%}
{% if 'application/vnd.vegalite.v2+json' in output.data %}
    <div id="vis{{cell['execution_count']}}"></div>
    <script type="text/javascript">
        var spec = {{ output.data['application/vnd.vegalite.v2+json'] }};
        var opt = {"renderer": "canvas", "actions": false};
        vegaEmbed("#vis{{cell['execution_count']}}", spec, opt);
    </script>
{% else %}
    {{super()}}
{% endif %}
{%- endblock data_priority -%}
