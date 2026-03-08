from(bucket: "{{ bucket }}")
  |> range(start: {{ start }}, stop: {{ stop }})
  |> filter(fn: (r) => r._measurement == "{{ measurement }}")
  |> filter(fn: (r) => r._field == "{{ field }}")
{% if extra_filters %}
  {% for field, value in extra_filters %}
  |> filter(fn: (r) => r["{{ field }}"] == "{{ value }}")
  {% endfor %}
{% endif %}
{% if aggregate %}
  |> {{ aggregate }}
{% endif %}