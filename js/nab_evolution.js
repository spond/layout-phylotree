//globals

var parse_date = d3.time.format("%Y%m%d");
var copy_re = /_([0-9]+)$/;

var show_date = function (d) { try {
                                if (date_annotations && date_annotations[d]) {return show_date_d (d) + " (" + date_annotations[d] + ")";};
                                return show_date_d (d);
                               }
                               catch (e) {
                                return d;
                               }
                             };
var show_date_axis = function (d) { try {
                                if (date_annotations && date_annotations[d]) {return date_annotations[d] + " (" + show_date_d(d) + ")";};
                                return show_date_d (d);
                               }
                               catch (e) {
                                return d;
                               }
                             };
                              
var show_date_d = d3.time.format("%b %Y");
var percentage = d3.format (".4g");
var percentage2 = d3.format (".2f");
var evolutionary_rates;
var _pos_data = [];
var _hxb2_map = [];
var _pos_sites = {};
var _pos_overall_data;
var _seq_data  = [];
var _ancestor_data = []; // one per date
var _mab_data  = [];
var _mab_features = [];
var _mrca_seq     = null;
var _seq_ids_to_dates = {};
var _all_sequences = {};
var _short_date_to_date = {};

var _tree_strings = {};
var _tree_widget = null;
var _tree_svg = null;

var _positive_selection = {};
var popover_obj = null;
var stored_clones = null;
var popover_feat_obj = null;
var stored_features = null;
var do_copy_number = false;
var date_annotations = null;

var _pos_sites_colors      = d3.scale.category10();
var _time_point_colors     = d3.scale.category10();
var _pos_sites_colors20      = d3.scale.category20();

var residues = "ACDEFGHIKLMNPQRSTVWY";
var _pos_sites_color_map = {};
for (var i = 0; i < residues.length; i++) {
    _pos_sites_color_map [residues[i]] = _pos_sites_colors20 (i); 
} 

_pos_sites_color_map ['ins'] = d3.rgb (128,128,128);
_pos_sites_color_map ['del'] = d3.rgb (128,128,128);

var _db_frequency_data_URL    = "/data/freqs.json";
var _db_frequency_information = null;
var _db_available_subtypes    = [];
var _toggle_evo_plot_checkers = 0;
var _positional_data_plot_height = 55;
var _hxb2_span = null;
var _alignment_length = null;

//handlers 

function get_copy_number (seq_name) {
    if (do_copy_number) {
        cpn = copy_re.exec(seq_name);
        if (cpn) {
            return parseFloat(cpn[1]);
        } else {
            return 1;
        }
    } else {
        return 1;
    }
}

function populateSelectWithValues (id, values, clear_first, secondary_values) {
	var select_element = d3.select ("#" + id);
	if (clear_first) {
		select_element.selectAll ("option").remove();
	}
	var option_set = select_element.selectAll ("option").data (values);
	
	option_set.enter ().append ("option");
	option_set.text (function (d) {return d;});
	if (secondary_values) {
	    option_set.attr ("value", function (d, i) {return secondary_values[i];});
	}
	option_set.exit ().remove();
			      
	return select_element;
}

function _tree_widget_node_colorizer (element, data) {
    element.style ("fill", _time_point_colors (_seq_ids_to_dates[data.name.toUpperCase()]));
}

function _tree_widget_show_date (data) {
    return _seq_ids_to_dates[data.name.toUpperCase()];
}

function set_handlers () {
    $('#_evo_plotme').button('loading');
    
    $('body').on('change', '#_pos_entropy', function (event) {
        event.stopPropagation(); // prevent default bootstrap behavior
        draw_positional_data ();
    });
    
    $('body').on('change', '#_subtype_for_conservation', function (event) {
        event.stopPropagation(); // prevent default bootstrap behavior
        draw_positional_data ();
    });

    $('body').on('change', '#_pos_entropy', function (event) {
        event.stopPropagation(); // prevent default bootstrap behavior
        draw_positional_data ();
    });

    $('body').on('change', '#_pos_conservation', function (event) {
        event.stopPropagation(); // prevent default bootstrap behavior
        draw_positional_data ();
    });

    $('body').on('change', '#_pos_sites_to_show', function (event) {
        event.stopPropagation(); // prevent default bootstrap behavior
        draw_positional_data ();
    });

    $('body').on('change', '#_site_subtype_for_conservation', function (event) {
        event.stopPropagation(); // prevent default bootstrap behavior
        site_id = $(this).data ('site');
        plot_site_conservation ([site_id,site_id], [270,300], "_seq_tertiary_plot", $(this).val())
        
    });
    
    $('#_seq_select_all').on('click', function (event) {
         event.stopPropagation(); // prevent default bootstrap behavior
        $("[data-seq-id]").prop ('checked', function (d, val) {
            return $(this).is(":visible") ? !val : val;
        }
        );      
    });
    
        

    
    $('#_tree_show_clone_names').on('change', function (event) {
    	_tree_widget.branch_name ($( this ).prop( "checked" ) ? null : _tree_widget_show_date ).update (true);
    });
    
    $('#_link_to_sequences').on('change', function (event) {
    	if (_tree_widget) {
    	    map_evolution_onto_tree ();
    	}
    });

    $('#_tree_regions_selector').on('change', function (event) {
    	var region = $(this).val();
    	
    	var option_values   = d3.keys (_tree_strings[region]).sort(),
    	    option_text = option_values.map (function (d) {var dp = parse_date.parse(d); if (dp) {return show_date (dp);} return d;});
    	
    	populateSelectWithValues ("_tree_time_point_selector", option_text, false, option_values);
    	$("#_tree_time_point_selector").trigger ("change");
    });

    $('#_tree_time_point_selector').on('change', function (event) {
    	var region = $("#_tree_regions_selector").val();
    	var date = $(this).val();
    	populateSelectWithValues ("_tree_scale_point_selector", d3.keys (_tree_strings[region][date]).sort());
    	$("#_tree_scale_point_selector").trigger ("change");
    });

	$('#_tree_scale_point_selector').on('change', function (event)  {
    	var region = $("#_tree_regions_selector").val();
    	var scale = $(this).val();
    	var date = $("#_tree_time_point_selector").val();
    	if (region && scale && date) {
    	    var tree_str = _tree_strings[region][date][scale];
    	    _tree_widget (tree_str).svg (_tree_svg).layout();
    	    map_evolution_onto_tree ();
    	}
	});

    $('body').on('click', '._seq_predef', function (event) {
        event.stopPropagation(); // prevent default bootstrap behavior
        var range = $(this).data ('range').split (',');
        $('#_seq_from').val(range[0]);
        $('#_seq_to').val(range[1]);
        $('#_seq_display').click();

    });

    $('body').on('click', '#_evo_plotme', function (event) {
        event.stopPropagation(); // prevent default bootstrap behavior
                
        plot_div_data (evolutionary_rates, $('#_evo_metrics').val(), 
                $("#_evo_metrics option:selected").get (). map (function (d, i) {return d.text;}),
                $('#_evo_regions').val(),"_evo_main_plot",[900,350],
                $("#_evo_metrics option:selected").get().map (function (d,i) {return [i,($(d).data("peer"))];}).filter (function (d) { return d[1] != undefined;})
                );


        plot_pheno_data (evolutionary_rates, $('#_evo_pheno').val(), 
                $("#_evo_pheno option:selected").get (). map (function (d, i) {return d.text;}),
                $('#_evo_regions').val(),"_evo_secondary_plot",[900,350]);
            
            
        return false;
    });

    $("#expand_spacing").on ("click", function (e) {
        _tree_widget.spacing_x (_tree_widget.spacing_x() + 1).update(true);
    });

    $("#compress_spacing").on ("click", function (e) {
        _tree_widget.spacing_x (_tree_widget.spacing_x() - 1).update(true);
    });

    function sort_nodes (asc) {
        _tree_widget.traverse_and_compute (function (n) {
                var d = 1;
                if (n.children && n.children.length) {
                    d += d3.max (n.children, function (d) { return d["count_depth"];});
                }
                n["count_depth"] = d;
            }); 
            _tree_widget.resort_children (function (a,b) {
                return (a["count_depth"] - b["count_depth"]) * (asc ? 1 : -1);
            });
    }

    $("#sort_original").on ("click", function (e) {
        _tree_widget.resort_children (function (a,b) {
            return a["original_child_order"] - b["original_child_order"];
        });
    });

    $("#sort_ascending").on ("click", function (e) {
        sort_nodes (true);
    });

    $("#sort_descending").on ("click", function (e) {
        sort_nodes (false);
    });
    
    function shift_coords (offset) {
        if (_hxb2_span) {
            $('#_seq_from').val(parseInt($('#_seq_from').val())+offset);
            $('#_seq_to').val(parseInt($('#_seq_to').val())+offset);
            $('#_seq_display').click();
        }
    }
    
    $("#_sequence_left1").on ("click", function (e) {
        e.stopPropagation();
        shift_coords (-1);
    });

    $("#_sequence_left5").on ("click", function (e) {
        e.stopPropagation();
        shift_coords (-5);
    });

    $("#_sequence_right1").on ("click", function (e) {
        e.stopPropagation();
        shift_coords (1);
    });

    $("#_sequence_right5").on ("click", function (e) {
        e.stopPropagation();
        shift_coords (5);
    });

    
    $('body').on('click', '#_seq_display', function (event) {
        event.stopPropagation(); // prevent default bootstrap behavior 
    
        function validate_range () {
            var range = [$('#_seq_from').val(),$('#_seq_to').val()].map (function (d) { return parseInt (d); });
            
            if (range[0] < 1) {
                range[0] = 1;
            }
            if (range[1] > _alignment_length) {
                range[1] = _alignment_length;
            }
            
            $('#_seq_from').val(range[0]);
            $('#_seq_to').val(range[1]);
            
            return range;
        }
    
        var seq_ids = {};
    
        $("[data-seq-id]").filter("input:checked").each (function (d) {
            seq_ids [$(this).attr ('data-seq-id')] = 1;
        }
        );
        
        _hxb2_span = validate_range().map (function (d) {return _hxb2_map.indexOf (d) + 1;});

        draw_sequence_data      (_hxb2_span,['_seq_positional_table','_seq_coordinates'], seq_ids); 
        if (_tree_widget) {
            map_evolution_onto_tree (_hxb2_span, seq_ids);
        }
        
        return false;
    });
    
}

function _update_evoplot_options () {

    var l = $('#_evo_regions').val().length;
    if (l > 1 && _toggle_evo_plot_checkers != 2) {
        $('#_evo_metrics').multiselect ('destroy');
        //d3.select ('#_evo_metrics').property ('multiple', false);
        $('#_evo_metrics').multiselect({'buttonText': function (options) { return _multiselect_button_label (options, 1, 1, 3);}});
        $('#_evo_pheno').multiselect ('destroy');
        //d3.select ('#_evo_pheno').property ('multiple', false);
        $('#_evo_pheno').multiselect(({'buttonText': function (options) { return _multiselect_button_label (options, 1, 1, 4);}}));
        _toggle_evo_plot_checkers = 2;
    }
    if (l <= 1 && _toggle_evo_plot_checkers != 1) {
        $('#_evo_metrics').multiselect ('destroy');
        //d3.select ('#_evo_metrics').property ('multiple', true);
        $('#_evo_metrics').multiselect({'buttonText': function (options) { return _multiselect_button_label (options, 1, 5, 3);}});
        $('#_evo_pheno').multiselect ('destroy');
        //d3.select ('#_evo_pheno').property ('multiple', true);
        $('#_evo_pheno').multiselect(({'buttonText': function (options) { return _multiselect_button_label (options, 1, 2, 4);}}));
        _toggle_evo_plot_checkers = 1;
    }
}

function _check_valid_evolplot () {
    var valid = false;
    try {
        var counts = [$('#_evo_regions').val().length,
                      $('#_evo_metrics').val().length,
                      $('#_evo_pheno').val().length];
        
        if (counts[0] > 1) {
            valid = counts[1] > 0 && counts[1] == 1 && counts[2] == 1;
            
        } else if (counts [0] == 1) {
            valid = counts[1] > 0 && counts[1] < 6 && counts[2] > 0 && counts[2] < 3;
        }
    }
    catch (e) {
        //console.log (e);
    }
    d3.select ("#_evo_plotme").attr ("disabled", valid ? null : true);
}

function _multiselect_button_label (options, min_required, max_allowed, limit, reg_exp, event_handler) {
    if (reg_exp == null) {
        reg_exp = /^(.+)$/;
    }
    
    text = '';
    
    if (options.length < min_required) {
        text = "<span class = 'alert alert-danger'>Select at least " + min_required + " value" + (min_required > 1 ? 's' : '') + "</span>";
    } else if (options.length > max_allowed) {
            text = "<span class = 'alert alert-danger'>Select no more than " + max_allowed + " value" + (max_allowed > 1 ? 's' : '') + "</span>";
        } else if (options.length == 0) {
              text = 'None selected ';
        }
        else if (options.length >= limit) {
          text = options.length + ' selected  ';
        } else {
        
            for (var i = 0; i < options.length; i++) {
                l_text = reg_exp.exec(options[i].label);
                if (l_text.length > 1) {
                    l_text = l_text [1];
                } else {
                    l_text = options[i].label;
                }
                text += l_text + ', ';
            }
            text = text.substr(0, text.length -2) + ' ';
        }
    
    if (event_handler != null) {
        event_handler ();
    } else {
        _check_valid_evolplot ();
    }

    return text + ' <b class="caret"></b>';
}

// loaders

function load_a_directory (name, dir, use_copy_numbers) {
    if (use_copy_numbers) {
        do_copy_number = true;
    }
    set_handlers();
    

 


/************** LOAD CASCADE ***********************/
    
    d3.json (dir + "dates.json", function (error, date_json) {
      
       if (date_json) {
           date_annotations = {};
           for (var d in date_json) {
            date_annotations [parse_date.parse (d)] = date_json[d];
           }
           _positional_data_plot_height += 25;
        }
       
       /************** LOAD RATES ***********************/
        d3.json (dir + 'rates.json', function (error, json) {
        for (k in json) {
                d = parse_date.parse(k);
                if (typeof json[k] == "string") {
                    rate_array = eval (json[k]);
                } else {
                    rate_array = json[k];
                }
                
                if (d) {
                    _pos_data.push ([d, rate_array]);
                    _positive_selection[d] = check_positive_selection(rate_array);
                } else {
                    _pos_overall_data = [k, rate_array];
                    _positive_selection['combined'] = check_positive_selection(rate_array);
                }
            }
        _pos_data.sort (function (a,b) {return a[0]-b[0];});   
            
        d3.json(dir + "frequencies.json", function(error, json) {
          if (error) return console.warn(error);
          
          _hxb2_mapper = [];
          
          for (var k in json) {
            _hxb2_mapper.push ([parseInt(k),parseInt(json[k]['HXB2'])]);
            if (get_site_residues(json,k).length > 1) {
                _pos_sites [+k] = json[k];
            }
          }
          _hxb2_mapper.sort (function (a,b) {return a[0]-b[0];});
          _hxb2_map = _hxb2_mapper.map (function (d) {return d[1];});
        });
       
       
        draw_positional_data ();
    });

       /************** LOAD SEQUENCES ***********************/
        d3.json (dir + 'sequences.json', function (error, json) {
    
       var flat_sequences    = [];
      
       for (k in json) {
            if (k == "MRCA") {
                _mrca_seq = json[k];
                flat_sequences.push (['mrca' , 'N/A', ""]);
                continue;
            }
            
            if (k == "Combined") {
                combined = [];
                for (i in json[k]) {
						combined.push ([i,json[k][i]]);
				}
            	_ancestor_data.push ([k,combined]);
            	continue;
            } else {
                date = parse_date.parse(k);
            }
            
            
            [['Observed', _seq_data], ['Ancestral', _ancestor_data]].forEach (function (d, ai) {
			        
					if (d[0] in json[k]) {
						sequence_source = json[k][d[0]];
			
						seq_for_date = [date, []];
						for (i in sequence_source) {
							if (ai == 0) {
							    var pd = parse_date.parse(k);
								_seq_ids_to_dates [i] = show_date(pd);
								_short_date_to_date[_seq_ids_to_dates [i]] = k;
								_all_sequences [i] = sequence_source[i];
								if (!_alignment_length) {
								    _alignment_length = _all_sequences [i].length;
								}
								
							}
							
							seq_for_date[1].push ([i,sequence_source[i]]);
						}
			
						d[1].push (seq_for_date);
					}
				}
			);
        }
        
                
        _seq_data.sort (function (a,b) {return a[0]-b[0];});
        _ancestor_data.sort (function (a,b) {return a[0]-b[0];});
 

       
        for (i = 0; i < _seq_data.length; i++) {
            d = _seq_data[i][0]
            sd = _seq_data [i][1];
            for (s = 0; s < sd.length; s++) {
                flat_sequences.push ([sd[s][0],show_date(d),""]);
            }
        }
               
        d3.select ("#_seq_select_table").selectAll ("tr")
          .data (flat_sequences)
          .enter ().append ("tr").selectAll ("td").data (function (d) {return d; })
          .enter ().append ("td").html (function (d, i, j) {if (i==2) {
            return '<input type="checkbox" data-seq-id = "' + flat_sequences[j][0] + '" checked/>';
          }
          return d;});
        
            $( '#_seq_select_limiter' ).on('input propertychange', function(event) {
            
                var filter_value = $(this).val().split (",");
                d3.select ("#_seq_select_table").selectAll ("tr").style ("display", function (d,i) { 
                    for (var k = 0; k < d.length; k++) {
                        var match_me = null;
                        if (typeof d[k] == "string") {
                            match_me = d[k]
                        } 
                        if (match_me) {
                            for (v = 0; v < filter_value.length; v++) {
                                if (match_me.indexOf(filter_value[v]) >= 0) {
                                    return "table-row";
                                }
                            }
                        }  
                    }
                    return "none";
                });
            }); 


    
       /************** LOAD TREES ***********************/
     	d3.json (dir + 'trees.json', function (error, json) {
            var tree_groups = 0;
            if (json) {
                for (var raw_date in json) {
                    date = (raw_date);
                    for (region in json[raw_date]) {
                        if (!(region in _tree_strings)) {
                            _tree_strings  [region] = {};
                        }
                        _tree_strings  [region] [date] = json[raw_date][region];
                        tree_groups++;
                    }
                }
            }
        

            if (tree_groups) {
                var width  = 800, //$(container_id).width(),
                    height = 600; //$(container_id).height()    
                
                _tree_widget = d3.layout.phylotree("body")
                    .size([height, width])
                    .separation (function (a,b) {return 0;})
                    .style_nodes (_tree_widget_node_colorizer)
                    .branch_name (_tree_widget_show_date);

                _tree_widget.node_span ('equal');
                
                if (do_copy_number) {
                     _tree_widget.node_span (function (a) { var m = copy_re.exec (a.name); try {return Math.sqrt(parseFloat (m[1]))} catch (e) {} return null;});
                    _tree_widget.options ({'draw-size-bubbles' : true}, false);               
                } else {
                    _tree_widget.options ({'draw-size-bubbles' : false}, false);
                }
                _tree_widget.options ({'selectable' : false}, false);


                _tree_svg = d3.select("#_tree_enclosure").append("svg")
                    .attr("width", width)
                    .attr("height", height);
                

            }
        
            populateSelectWithValues ("_tree_regions_selector", d3.keys (_tree_strings).sort());
            $("#_tree_regions_selector").val("gp160").trigger ("change");
        }
        );

    
        d3.json (dir + 'mab.json', function (error, json) {
        
        	if (!json) {
        		return;
        	}
        	       
            var mab_names         = [],
                sequences_by_date = {},
                sequences_by_date_array = [],
                unique_dates      = {};
                
            _mab_data = [];
            
            for (var mname in json) {
                mab_names.push (mname);
                json[mname]['predictions'].forEach (function (seq) {
                    var seq_id = seq['id'];
                    var d = _seq_ids_to_dates[seq_id];
                    if (d) {
                        d = _short_date_to_date[d];
                        if (! (d in sequences_by_date)) {
                            sequences_by_date[d] = {};
                        }
                    
                        if (! (d in unique_dates)) {
                            unique_dates [d] = 1;
                        }
                    
                        if (!(mname in sequences_by_date[d])) {
                            sequences_by_date[d][mname] = [];
                        }
                    
                        sequences_by_date[d][mname].push([seq_id, seq['value'] > 0 ? true : false,seq['features']]);
                    }
                }
                );
            }
                        
         
            sequences_by_date_array.sort (function (a,b) {return a[0]-b[0];});
            
            
            var unique_dates_d = [];
            
            for (var kk in unique_dates) {
                unique_dates_d.push(kk);
            }
            
            unique_dates_d.sort();
            
            for (var i = 0; i < unique_dates_d.length; i ++) {
                var date = unique_dates_d[i];
                sequences_by_date_array.push ([date, sequences_by_date[date]])
            }

             
            sequences_by_date_array = sequences_by_date_array.map (function (d,i) { 
                var string_date = d[0];
                var mab_array = [];
                mab_names.forEach (function (mn) {
                        mab_array.push (d[1][mn]);
                    }
                );
                return [string_date, mab_array];
            });
            
            
            mab_names.forEach (function (mname, i) {
                var entry = [mname];
                sequences_by_date_array.forEach (function (d) {
                        entry.push (d[1][i]);
                    }
                ); 
                _mab_data.push (entry);
            });
            
            
            _mab_data.sort ((function (a,b) {return a[0]-b[0];}));
            unique_dates_d.splice (0, 0, "mAb");
            d3.select ('#_nab_table_head').selectAll ("td").data (unique_dates_d).enter().append ("td")
                .text (function (d,i) {return i > 0 ? show_date(parse_date.parse(d)) : d;});
                
            var red_white = d3.interpolateRgb("white", "red");   
                    
            var feature_lists = [];        
                
            function bin_features (list) {
                var feat_counts = {},
                    feat_array = [];
                    
                //console.log (list);
                    
                list.forEach (function (d) {
                    
                    d.forEach (function (l) {
                        if (l in feat_counts) {
                            feat_counts [l] ++;
                        } else {
                            feat_counts [l] = 1;
                        }
                    });
                });
                
 
                for (f in feat_counts) {
                    feat_array.push ([f,feat_counts[f]]);
                }    
                
                //feat_array.sort ((function (a,b) {return a[0]-b[0];}));
                return feat_array;
            }    
                
            _mab_features = [];    
                
            d3.select ('#_nab_table_body').selectAll ("tr").data (_mab_data).enter().append ("tr")
                .selectAll ('td').data (function (d) {return d;}).enter().append ("td")
                .html (function (d,i,j) { 
                    if (i == 0) return d; 
                    
                    if (_mab_features.length == j) {
                        _mab_features.push([]);
                    }
                    
                    
                    var resistant   = d.filter (function (e,i) { return e[1];}),
                        susceptible = d.filter (function (e,i) { return !e[1];});
                    
                    _mab_features[j].push ([bin_features (susceptible.map (function (d) {return d[2];})),
                                          bin_features (resistant.map (function (d) {return d[2];}))]);
                    
                    this.frac = resistant.length/d.length;
                    return "<a href = '#' class = '_nab_cells' data-cell-id = '" + j + "-" + (i-1) + "'>" 
                            + percentage(100*this.frac) + "%</a>";
            })
                .style ("background-color", 
                    function (d) {
                        return red_white (this.frac);
                    });                    
                            
            $('._nab_cells').on('click', function (event) {
                display_feature_info ($(this), $(this).attr ("data-cell-id").split('-'));
               
                //console.log (_mab_features[parseInt(idx[0])][parseInt (idx[1])]);
            });
            
        });   



    });
    });


/************** LOAD BASELINE FREQS ***********************/

    d3.json (_db_frequency_data_URL, function (error, json) {
        if (json) {
            if ("overall" in json && "subtype" in json) {
                _db_frequency_information = json;
                for (subtype in json["subtype"]) {
                    subtype_counts = json["subtype"][subtype];
                    if (_add_dict_values(subtype_counts[0]) > 100) {
                        _db_available_subtypes.push (subtype);
                    }
                }
                _db_available_subtypes = _db_available_subtypes.sort();
                _db_available_subtypes.splice (0, 0, 'Combined');
                
                add_options_to_select ("_subtype_for_conservation", _db_available_subtypes.map (function (d) {return [d,d];}));
                add_options_to_select ("_site_subtype_for_conservation", _db_available_subtypes.map (function (d) {return [d,d];}));
                

            }
        }
        //console.log (error);
    });

/************** LOAD RATES/PHENO ***********************/

    d3.tsv(dir + "rates_pheno.tsv")
    .row(function(d) { return d; })
    .get(function(error, rows) { 
        regions = {};
        rows.forEach (function (d) {
            d.Date = parse_date.parse(d.Date);
            regions[d.Segment] = 1;
            for (k in d) {
                if (k != "Segment" && k != "Date") {
                    d[k] = +d[k];
                }
            }
        });
        
        
        d3.select ("#_evo_regions").selectAll ("option").data(d3.keys (regions))
        .enter()
        .append ("option")
        .text (function (d) {return d;});
        
        d3.select ("#_evo_regions").selectAll ("option").property ('selected', function (d,i) {return i == 0;});
        
        

        
        var pheno = [['Length', 'Sequence Length'], 
                     ['PNGS', "Number of PNGS"],
                     ["IsoelectricPoint", "Isoelectric Point"]];
                     
        var geno  = [   
                        ['ds_divergence', 'dS divergence'],
                        ['dn_divergence', 'dN divergence'],
                        ['total_divergence', 'Divergence'],
                        ['ds_diversity', 'dS diversity', 'ds_divergence'],
                        ['dn_diversity', 'dN diversity', 'dn_divergence'],
                        ['total_diversity', 'Diversity', 'total_divergence']
                    ];       
                    
        d3.select ("#_evo_pheno").selectAll ("option").data(pheno)
         .enter()
         .append ("option")
         .attr ("value", function (d) { return d[0]; })
         .text (function (d) {return d[1];});        
    
        d3.select ("#_evo_metrics").selectAll ("option").data(geno)
         .enter()
         .append ("option")
         .attr ("value", function (d) { return d[0]; })
         .text (function (d) {return d[1];})
         .attr ("data-peer", function (d) { if (d.length < 3 || d[2].length == 0) return null; return d[2]});
 
 
        d3.select ("#_evo_regions").selectAll ("option").property ('selected', function (d,i) {return i == 0;});
        //d3.select ("#_evo_regions").selectAll ("option").property ('selected', function (d,i) {return i == 1 || i == 3;});
        $('#_evo_regions').multiselect({'buttonText': function (options) { return _multiselect_button_label (options, 1, 5, 4);}, 'onChange': _update_evoplot_options });

        //d3.select ("#_evo_metrics").selectAll ("option").property ('selected', function (d,i) {if (i == 2) return true; return false;});
        
        d3.select ("#_evo_metrics").selectAll ("option").property ('selected', function (d,i) {if (i < 3 || i == 5) return true; return false;});
        $('#_evo_metrics').multiselect({'buttonText': function (options) { return _multiselect_button_label (options, 1, 5, 3);}});
        d3.select ("#_evo_pheno").selectAll ("option").property ('selected', function (d,i) {if (i < 2) return true; return false;});
        //d3.select ("#_evo_pheno").selectAll ("option").property ('selected', function (d,i) {if (i < 1) return true; return false;});
        $('#_evo_pheno').multiselect(({'buttonText': function (options) { return _multiselect_button_label (options, 1, 2, 4);}}));
       
        _check_valid_evolplot();
        $('#_evo_plotme').button('complete');
        evolutionary_rates = rows;
        $('#_evo_plotme').click();
       
        d3.select ("#allofit").style ('display', 'block');      
        d3.select ("#loading_bar").style ('display', 'none');      
     });       
    
}

var global_summary_filter = function (d,i) {return true;};

function map_evolution_onto_tree_get_range (name, span, leaves, inodes) {
    if (name == 'root') {
        name = 'Node0';
    }
    if (name in leaves) {
        return leaves[name].slice (span[0]-1, span[1]);
    } 
    if (name in inodes) {
        return inodes[name].slice (span[0]-1, span[1]);
    }
    return null;
}

function array2dict (a,d) {
    for (i = 0; i < a.length; i++) {
        d[a[i][0]] = a[i][1];
    }
}

function map_evolution_onto_tree (span, visible_clone_list) {

    span = span || _hxb2_span;
    
    var tree_date = parse_date.parse($("#_tree_time_point_selector").val());
    
    if (tree_date) {    
        for (k = 0; k < _seq_data.length; k++) {
            if (_seq_data[k][0] - tree_date == 0) {
                break;
            }
        }
    } else {
        tree_date = $("#_tree_time_point_selector").val();
        for (k = 0; k < _ancestor_data.length; k++) {
            if (_ancestor_data[k][0] == tree_date) {
                break;
            }
        }
    }
    
    var inodes = {};
    
    array2dict ( _ancestor_data [k][1], inodes);
              
    var show_subs = $( "#_link_to_sequences" ).prop( "checked" ) && span;
    
     
    _tree_widget.style_edges (
        function (element, data) {
            $(element[0][0]).tooltip('destroy');
            
            if  (show_subs) {
                var source_seq = map_evolution_onto_tree_get_range(data.source.name, span, _all_sequences, inodes),
                    target_seq = map_evolution_onto_tree_get_range(data.target.name, span, _all_sequences, inodes);
            
                if (source_seq && target_seq) {
                    var diffs = 0;
                    var html = "";
                    for (k = 0; k < source_seq.length; k++) {  
                        if (source_seq[k] != target_seq[k]) {
                            if (diffs) {
                                html += "<br/>";
                            }
                            //console.log (span);
                            html += "" + _hxb2_map[span[0] + k -1] + " " + source_seq[k] + "&rarr;" + target_seq[k];
                            diffs +=  1;
                        }
                    }
                
                
                    if (diffs) {
                        $(element[0][0]).tooltip({'title' : html, 'html' : true, 'trigger' : 'hover', 'container' : 'body', 'placement' : 'auto'}); 
                    } 
                    element.style ("stroke-width", Math.min (1+diffs,10)).style ("stroke", diffs>0 ? "red" : null);
                } else {
                    //console.log (data);
                }
            } else {
                element.style ("stroke-width", 1).style ("stroke", null);
            }
        });
        
    _tree_widget.update (true);
}

function draw_sequence_data (span, id, visible_clone_list) {  

    var show_selected = $( "#_seq_selected" ).prop( "checked" );
    
    
    var pngs_re = /N\-*[^P]\-*[ST]\-*[^P]/g;
    
    function display_sequence (seq_info,r,c) {
        new_s = "";
        
        
        pngs_indices = [];
        while ( (result = pngs_re.exec(seq_info[0])) ) {
            pngs_indices.push([result.index, result.index + result[0].length-1]);
        }
        
        current_pngs = 0;
        
        for (i = 0; i < seq_info[0].length; i++) {
            if (current_pngs < pngs_indices.length) {
                if (i == pngs_indices[current_pngs][0]) {
                    new_s += "<span class = 'pngs'>";
                } 
            }
            if ((r > 0 || c > 0) && seq_info [0][i] == baseline_sequence[i]) {
                new_s += ".";
            } else {
                new_s += "<span style = 'color: " +_pos_sites_color_map[seq_info[0][i]] + "'>" + (seq_info[0][i] == " " ? "&nbsp;" : seq_info[0][i])  + "</span>";
            }
            if (current_pngs < pngs_indices.length) {
                if (i == pngs_indices[current_pngs][1]) {
                    new_s += "</span>";
                    current_pngs++;
                }
                
            }
        }
        return new_s + (seq_info [1][0] > 0 ? ("&nbsp;<a href = '#' class = 'clone_info' data-sequence-list = '" + seq_info[1][1].join(';') + "'>" + seq_info[1][0] + "</a></span> (" + percentage(100*seq_info[1][0]/seq_info[2]) + "%)") : "");
    }
        
    var sorted_sequences_by_date = [];
    var baseline_sequence = null;
    
    var show_mrca = 'mrca' in visible_clone_list;
    
    if (show_mrca && _mrca_seq != null) {
        baseline_sequence = _mrca_seq.slice(span[0]-1,span[1]);
    }
    
    var total = 0,
        shown = 0;
   
    for (i = 0; i < _seq_data.length; i++) {
        var unique_strings = {};
        var local_total = 0;
        for (s = 0; s < _seq_data[i][1].length; s++) {
            cpn = get_copy_number (_seq_data[i][1][s][0]);
            local_total += cpn;
            if (visible_clone_list && !(_seq_data[i][1][s][0] in visible_clone_list)) {
                continue;
            } 
            shown += cpn;
            ss = _seq_data[i][1][s][1].slice(span[0]-1,span[1]);
            if (ss in unique_strings) {
                unique_strings[ss][0] += cpn ;
                unique_strings[ss][1].push (_seq_data[i][1][s][0]);
            } else {
                unique_strings[ss] = [cpn,[_seq_data[i][1][s][0]]];
            }
        }
        total += local_total;
        
        
        sorted_unique_strings = []
        for (s in unique_strings) {
            sorted_unique_strings.push ([s, unique_strings[s],local_total]);
       }
     
        if (sorted_unique_strings.length) {
        
            sorted_unique_strings.sort (function (a,b) {return b[1][0]-a[1][0];});
            if (show_selected) {
                selected_string = "";
                for (var site = span[0]; site <= span[1]; site++) {
                    selected_string += (_positive_selection[_seq_data[i][0]].indexOf (site)>=0) ? "+" : " ";
                }
                if (selected_string.indexOf('+') >= 0) {
                    sorted_unique_strings.push ([selected_string,[-1,0],0]);
                }
            }
        
            if (baseline_sequence == null && sorted_sequences_by_date.length == 0) {
                baseline_sequence = sorted_unique_strings[0][0];
            }
            sorted_sequences_by_date.push ([show_date(_seq_data[i][0]), sorted_unique_strings]);
        }
        
    }
    
    if (show_mrca && _mrca_seq != null) {
        sorted_sequences_by_date.splice (0,0,['MRCA', [[baseline_sequence, [1, ['mrca']],1]]]);
    }
        
    //console.log (_seq_data, sorted_sequences_by_date);
    
    d3.select ('#_seq_selected_count').text (shown);
    
    d3.select ('#' + id[1]).html (function () {
        coordinates = ["","",""];
        selected_string = "";
        var last_hs = -5;
        for (var s = span[0]; s <= span[1]; s++) {
            var hs = _hxb2_map[s-1];
            if (hs < 10) {
                str = "  " + hs;
            } else {
                if (hs < 100) {
                    str = " " + hs;
                } else {
                    str = "" + hs;
                }
           }
           if (show_selected) {
                selected_string += (_positive_selection['combined'].indexOf (s)>=0) ? "+" : "&nbsp;";
           }  
           
           if (last_hs == hs) {
            str = "INS";
           }
                 
           coordinates[0] += "<span class = '_seq_hover_seq' data-coord = '" + s +"'>" + (str[0] == " " ? "&nbsp;" : str[0]) + "</span>";
           coordinates[1] += "<span class = '_seq_hover_seq' data-coord = '" + s +"'>" + (str[1] == " " ? "&nbsp;" : str[1]) + "</span>";
           coordinates[2] += "<span class = '_seq_hover_seq' data-coord = '" + s +"'>" + (str[2] == " " ? "&nbsp;" : str[2]) + "</span>";

           last_hs = hs;
        }
        return coordinates[0] + "<br/>" + coordinates[1] + "<br/>" + coordinates[2] + (show_selected ? "<br/>" + selected_string : "");
    });
    
    $('._seq_hover_seq').mouseenter(function (event) {
        event.stopPropagation();
        site_id = parseInt($(this).data ('coord'));
        plot_site_data ([site_id,site_id], [270,200], "_seq_secondary_plot");
        $("#_site_subtype_for_conservation").data ('site', site_id);
        $("#_site_subtype_conservation_div").show();
    
        plot_site_conservation ([site_id,site_id], [270,300], "_seq_tertiary_plot", $("#_site_subtype_for_conservation").val());
        
    });        
    var seq_table = d3.select ('#' + id[0]);

    
    seq_table.selectAll ('tr').remove();
    seq_table.selectAll ('tr').data (sorted_sequences_by_date).enter()
                  .append ('tr').selectAll ('td').data (function (d) {return d; }).enter()
                  .append ('td').text (function (d,i) {if (i == 0) return d; return ""; }).
                  attr ('class',function (d, i) { if (i==1) return "sequence"; else return "seq-date"});

    seq_table.selectAll ('.sequence').selectAll ('ul').
            data (sorted_sequences_by_date.map (function (d) {return d[1];})).
            enter().append ('ul').attr ('class', 'list-unstyled').
            data (function (d) {return [d];}).selectAll('li').
            data (function (d) {return d;}).enter().append ('li').
            html (display_sequence);
            
    $( '.clone_info' ).on( 'click', function(event) {
        display_clone_info ($(this), $(this).attr ("data-sequence-list"));
    });
            
}

function generate_clone_info_html (clone_list) {
    var clones = clone_list.split (";");
    clones.sort();
    var html = "<table class = 'table-striped table table-hover'><thead><tr><th>Sequence ID</th></tr><tbody>";    
    for (k = 0; k < clones.length; k++) {
        html += "<tr><td>" + clones[k] + "</td></tr>";
    }
    return html + "</tbody></table>";
}


function generate_feature_info_html (feat_list) {
 
    var html = "<table class = 'table-striped table table-hover'><thead><tr><th>Feature</th><th>Count (S/R)</th></tr><tbody>";    
    
    var combined_features = {};
    
    feat_list.forEach (function (s,i) {
        s.forEach (function (feat) {
            if (!(feat[0] in combined_features)) {
                combined_features[feat[0]] = [0,0];
            }  
            combined_features[feat[0]][i] += feat[1];
        });
        
    });
    
    for (k in combined_features) {
       html += "<tr><td>" + k + "</td><td>" + combined_features[k][0] + "/" + combined_features[k][1] + "</td></tr>";
    }
    return html + "</tbody></table>";
}

function display_clone_info (obj, selected_clones) {
   if (popover_obj) {
        popover_obj.popover ('destroy');
   }
   if (stored_clones == selected_clones) {
    popover_obj = null;
    stored_clones = null;
    return;
   }
   obj.popover({animation : false, placement: "left", title: null, html: true, content: generate_clone_info_html(selected_clones), trigger : "click"});
   obj.popover('show');
   popover_obj = obj;
   stored_clones = selected_clones;
}

function display_feature_info (obj, selected_feats) {
   if (popover_feat_obj) {
        popover_feat_obj.popover ('destroy');
   }
   if (stored_features == selected_feats) {
        popover_feat_obj = null;
        stored_features = null;
        return;
   }
   var idx = selected_feats.map (function (d) {return parseInt(d);});
   obj.popover({animation : false, placement: "left", title: null, html: true, 
            content: generate_feature_info_html(_mab_features[idx[0]][idx[1]]), trigger : "click"});
   obj.popover('show');
   popover_feat_obj = obj;
   stored_features = selected_feats;
}


function handle_positional_data (overall, positional, id, dim, labels, annotator) {
    var two_d = overall[1][0].length > 1;
    var show_selected = $( "#_pos_selected" ).prop( "checked" );
    
    function brushed() {
      x.domain(brush.empty() ? x_overall.domain() : brush.extent());
      for (k = 0; k < focus_plots.length; k++) {
          focus_plots[k].select("._pos_dS").attr("d", area_objects[k][0]);
          if (two_d) {
            focus_plots[k].select("._pos_dN").attr("d", area_objects[k][1]);
          }
          focus_plots[k].selectAll (".selected_site").
                attr ("cx", function (d) {return x(d);});
                
      }
      svg.select(".pos.x.axis").call(xAxis);
      if ($( "#_pos_conservation" ).prop( "checked" )) {
        plot_site_conservation (x.domain(), [300,200], "_pos_secondary_plot", $('#_subtype_for_conservation').val());    
      }    
      else {
        plot_site_data (x.domain(), [300,150], "_pos_secondary_plot");
      }
    }
    
    var max_values = [d3.max (overall[1].map (function (d) { return d[0]; })),two_d ? d3.max (overall[1].map (function (d) { return d[1]; })) : 0];
    
    for (k = 0; k < positional.length; k += 1) {
        mx = d3.max (positional[k][1].map (function (d) { return d[0]; }));
        max_values[0] = Math.max (max_values[0], mx);
        if (two_d) {
            mx = d3.max (positional[k][1].map (function (d) { return d[1]; }));
            max_values[1] = Math.max (max_values[1], mx);
        }
    }
    
        
    var margin              = {top: 20, right: 0, bottom: 50, left: 30},
        width               = dim[0] - margin.left - margin.right,
        height              = dim[1]*(positional.length+1);

    var svg = d3.select("#" + id)
        .attr("width", width + margin.left + margin.right)
        .attr("height", height + margin.top + margin.bottom);
    
    svg.selectAll("path").remove();
    svg.selectAll("g").remove();   
    svg.selectAll("defs").remove();   
    
    svg.append("defs").append("clipPath")
        .attr("id", "clip")
        .append("rect")
        .attr("width", width)
        .attr("height", height); 
        
    svg = svg.append ("g");
    

    var x = d3.scale.linear()
        .range([0, width]).domain ([1,positional[0][1].length]);
 
    var x_overall = d3.scale.linear()
        .range([0, width]).domain ([1,positional[0][1].length]);
        
    var brush = d3.svg.brush()
        .x(x_overall)
        .on("brush", brushed);    

    var y = d3.scale.linear()
        .range([dim[1], 0]).domain ([two_d ? -max_values[1] : 0,max_values[0]]);
    
    var xAxis = d3.svg.axis()
        .scale(x)
        .orient("bottom");

    var xAxis_overall = d3.svg.axis()
        .scale(x_overall)
        .orient("top");

    
    var focus_plots  = [];  
    var area_objects = [];
      
    for (plot_id = 0; plot_id <= positional.length; plot_id++) {
    
        plot_svg = svg.append ("g")
                   .attr("transform", "translate(" + margin.left + "," + (margin.top + dim[1]*plot_id) + ")");
    
        
        var the_data = plot_id > 0 ? positional[plot_id-1] : overall;
        
        var local_areas = [];
        
        local_areas[0] = d3.svg.area()
            .x(function(d,i) { return x(i+1); })
            .y0(function(d) {  return y(0); })
            .y1(function(d) {  return y(d[0]); })
            .interpolate('step');
            
        if (two_d) {    
            local_areas[1] = d3.svg.area()
                .x(function(d,i) { return x(i+1); })
                .y0(function(d) {  return y(-d[1]); })
                .y1(function(d) {  return y(0); })
                .interpolate('step');
        }
            
        area_objects.push (local_areas);
            
        plot_svg.append ("text")
              .attr("transform", "translate(0," + y(0) + ") rotate (-90)")
              .attr("y", "0")
              .attr("dy", "-0.5em")
              .style("text-anchor", two_d ? "middle" : "start")
              .text(plot_id ? show_date(the_data[0]) : the_data[0]);

        plot_svg.append("path")
              .datum(the_data[1])
              .attr("class", "_pos_dS")
              .attr("clip-path", "url(#clip)")
              .attr("d", local_areas[0]);

        if (two_d) {
            plot_svg.append("path")
                  .datum(the_data[1])
                  .attr("clip-path", "url(#clip)")
                  .attr("class", "_pos_dN")
                  .attr("d", local_areas[1]);
        }
              
        if (plot_id == 0) {
            plot_svg.append("g")
              .attr("class", "x brush")
              .call(brush)
            .selectAll("rect")
              .attr("y", -6)
              .attr("height", dim[1] + 7);
              
                
        } else {
            focus_plots.push (plot_svg);

        }
        
        if (show_selected) {
            plot_svg.selectAll (".selected_site").
                data (_positive_selection[plot_id>0 ? the_data[0] : "combined"]).
                enter().
                append("circle").
                attr ("cx", function (d) {return x(d);}).
                attr ("cy", function (d) {return y(0);}).
                attr ("r", 4).
                attr ("class", "selected_site").
                append("title").
                text (function (d) { return "Codon " + d; });

        }
    }

    svg.append("g")
      .attr("class", "pos x_overall axis")
      .attr("transform", "translate("+ margin.left + "," + (margin.top-2) + ")")
      .call(xAxis_overall)
      .append("text")
      .attr("transform", "translate(0,0)")
      .attr("dy", "-.4em")
      .attr("dx", ".25em")
      .style("text-anchor", "start")
      .text("Site");    

    svg.append("g")
      .attr("class", "pos x axis")
      .attr("transform", "translate("+ margin.left + "," + (margin.top+height) + ")")
      .call(xAxis)
      .append("text")
      .attr("transform", "translate(0,0)")
      .attr("dy", "+1.1em")
      .attr("dx", ".25em")
      .style("text-anchor", "start")
      .text("Site");    

     if (labels) {
     
        var legend_dim = {x: 20, y:height + margin.top + margin.bottom/3, spacer:15, margin:5, font: 10, x_step : 100};
        
        var me_colors = ['#2E66FF', '#FFB314'];
    
        var legend = svg.append("g")
          .attr("class", "_evo_legend")
          .attr("x", legend_dim.x)
          .attr("y", legend_dim.y)
          .attr("transform", "translate("+legend_dim.x+","+legend_dim.y+")");
        
        legend.selectAll('g').data(labels)
          .enter()
          .append('g')
          .each(function(d, i) {
            var g = d3.select(this);
            g.append("rect")
              .attr("y", legend_dim.spacer)
              .attr("x", i*(legend_dim.spacer + legend_dim.margin + legend_dim.x_step))
              .attr("width", legend_dim.spacer)
              .attr("height", legend_dim.spacer)
              .style("fill", me_colors[i]);
        
            g.append("text")
              .attr("y", 2*legend_dim.spacer)
              .attr("x", i*(legend_dim.spacer + legend_dim.margin + legend_dim.x_step) + legend_dim.spacer + legend_dim.font/4)
              .style("fill", me_colors[i])
              .text(function (d) {return d;});
          
          });
     }

      if ($( "#_pos_conservation" ).prop( "checked" )) {
        plot_site_conservation (x_overall.domain(), [300,200], "_pos_secondary_plot", $('#_subtype_for_conservation').val());
      }    
      else {
        plot_site_data (x_overall.domain(), [300,150], "_pos_secondary_plot");
      }
}
 
 
//plotters

function plot_site_conservation (site_range, dim, id, subtype) {


    var _max_sites_to_plot = parseInt($( "#_pos_sites_to_show" ).val());

    var margin = {top: 5, right: 30, bottom: 10, left: 30},
        buffer = 10,
        width  = dim[0] - margin.left - margin.right,
        height = dim[1] - 2*buffer;


    var svg = d3.select("#" + id)
        .attr("width", width + margin.left + margin.right)
        .attr("height", (dim[1])*_max_sites_to_plot + margin.top + margin.bottom);

    svg.selectAll("path").remove();
    svg.selectAll("g").remove();    
    
    var sites_plotted = 0;
   

    for (site_id = Math.ceil(site_range[0]); site_id <= Math.floor(site_range[1]) && sites_plotted < _max_sites_to_plot; site_id++) {
        
        var residue_counts = get_db_residue_counts (_hxb2_map [site_id-1] - 2, subtype);
        
        if (residue_counts) {
            var total = _add_dict_values (residue_counts, {'ins' : 1});
            var conservation_data = [];
        
        
            for (k in residue_counts) {
                conservation_data.push ([k, residue_counts[k] / total * 100.]);
            }
        
            conservation_data = conservation_data.sort (function (a,b) {return b[1]-a[1];});
        

            var plot_svg = svg.append("g")
                .attr("transform", "translate(" + margin.left + "," + (margin.top + (dim[1] + buffer) * sites_plotted)  + ")");
       
            sites_plotted += 1;
 
            var x = d3.scale.linear().domain ([1e-5,100])
                .range([0, width]);

            var y = d3.scale.ordinal()
                .rangeRoundBands([height, 0],0.1).domain (conservation_data.map (function (d) {return d[0];}));
            
        

            var xAxis = d3.svg.axis()
                .scale(x)
                .orient("bottom");
                //.ticks(10, function(d) { ex = Math.round (Math.log (d) / Math.LN10); if (ex >= 0) {return d.toFixed (0);} if (ex >= -4) {return d.toFixed (-ex);}; return d.toExponential (1); });
    
            var yAxis = d3.svg.axis()
                .scale(y)
                .orient("left");
           
            plot_svg.selectAll(".bar")
                  .data(conservation_data)
                  .enter().append("rect")
                  .attr("class", "bar")
                  .attr("y", function(d) { return y(d[0]); })
                  .attr("height", y.rangeBand())
                  .attr("x", function(d) { return 0; })
                  .attr("width", function(d) { return x(d[1]); });

            //console.log (y.rangeBand());

            plot_svg.selectAll(".barlabel")
                  .data(conservation_data)
                  .enter().append("text")
                  .attr("class", "barlabel")
                  .attr("y", function(d) { return y(d[0]) + 0.5*y.rangeBand (); })
                  .attr("dy", "0.33em")
                  .attr("x", function(d) { return x(d[1]) + 2; })
                  .text(function(d) { return percentage2(d[1]); }); 
              
            plot_svg.append ("text")
                .attr ("class", "barlabel")
                .attr ("x", width)
                .attr ("y", 15)
                .style("text-anchor", "end")
                .text ("HXB2 site " + (_hxb2_map [site_id-1] - 1));  
              
            plot_svg.append("g")
              .attr("class", "x axis")
              .attr("transform", "translate(0," + height + ")")
              .call(xAxis)
                  .append ("text")
                  .attr("transform", "translate(" + width + "," + 0 + ")")
                  .attr("dy", "3em")
                  .attr("dx", "-.4em")
                  .style("text-anchor", "end")
                  .attr ("class","axis")
                  .text ("Frequency, %");

              plot_svg.append("g")
                  .attr("class", "y axis")
                  .call(yAxis);
        }
              
    }

    svg.attr("height", (dim[1]+buffer)*sites_plotted + margin.top + margin.bottom);
}


function plot_site_data (site_range, dim, id) {

    _max_sites_to_plot = parseInt($( "#_pos_sites_to_show" ).val());
    
    var margin = {top: 80, right: 50, bottom: 10, left: 30},
        buffer = 5,
        width  = dim[0] - margin.left - margin.right,
        height = dim[1] - 2*buffer;

    var x = d3.time.scale()
        .range([0, width]);

    var y = d3.scale.linear()
        .range([height, 0]);

    var xAxis = d3.svg.axis()
        .scale(x)
        .orient("top");

    var yAxis = d3.svg.axis()
        .scale(y)
        .orient("left");

    var svg = d3.select("#" + id)
        .attr("width", width + margin.left + margin.right)
        .attr("height", (dim[1])*_max_sites_to_plot + margin.top + margin.bottom);

    svg.selectAll("path").remove();
    svg.selectAll("g").remove();    

    frequency_data   = [];
    
    var sites_plotted = 0;
    var all_residues = {};

    for (site_id = Math.ceil(site_range[0]); site_id <= Math.floor(site_range[1]) && sites_plotted < _max_sites_to_plot; site_id++) {
        if (site_id in _pos_sites) {
            var stack_data   = [];
    
            my_residues = get_site_residues (_pos_sites, site_id);

            for (var r in my_residues) {
                stack_data.push ({'residue' : my_residues[r], 'values' : []});
                all_residues [my_residues[r]] = 1;
            }
        
            var dates = [];
            var sums  = [];
        
            for (var k in _pos_sites[site_id]) {
                if (k != 'HXB2') {
                    var md = parse_date.parse(k);
                    dates[md] = 1;
                    var sum = 0;
                    for (var r in my_residues) {
                        if (my_residues[r] in _pos_sites[site_id][k]) {
                            var count = _pos_sites[site_id][k][my_residues[r]];
                            stack_data[r]['values'].push ({'x' : md, 'y' : count});
                            sum += count;           
                        } else {
                            stack_data[r]['values'].push ({'x' : md, 'y' : 0});
                        }
                    }
                    sums.push (sum);
                    dates.push (md);
                }
            }        
    
            x.domain (d3.extent (dates));
            xAxis.tickValues (dates).tickFormat (show_date_axis);
            y.domain ([0,d3.max (sums)]);
    
    
            var stack = d3.layout.stack()
            .values(function(d) { return d.values; });
     
            var area = d3.svg.area()
            .x(function(d) { return x(d.x); })
            .y0(function(d) { return y(d.y0); })
            .y1(function(d) { return y(d.y0 + d.y); })
            .interpolate ('monotone');      
    
    
            var plot_svg = svg.append("g")
            .attr("transform", "translate(" + margin.left + "," + (margin.top + dim[1] * sites_plotted)  + ")");

            
            plot_svg.selectAll("path")
                .data(stack(stack_data))
                .enter().append("path")
                .attr("class", "area")
                .style("fill", function (d) { return _pos_sites_color_map[d.residue]; })
                .attr("d", function(d) {  return area(d.values); })
                .append ("title").text (function (d) {return d.residue});
            
            plot_svg.append ("text")
              .attr("transform", "translate(" + width + "," + height + ")")
              .attr("dy", "-.4em")
              .attr("dx", "-.4em")
              .style("text-anchor", "end")
              .attr ("class","axis")
              .text("" + site_id + " (HXB2 " + _pos_sites[site_id]['HXB2'] + ")");

            if (sites_plotted == 0) {
                svg.append("g")
                  .attr("class", "x axis")
                  .attr("transform", "translate(" + margin.left + "," + margin.top + ")")
                  
                  .call(xAxis).selectAll("text").attr ("transform", "rotate(-90)").attr("dx","0.5em").attr ("dy", "1.25em")
                              .style("text-anchor", "start");


            }     
            
            plot_svg.append("g")
              .attr("class", "y axis")
              .attr("transform", "translate(0,0)")
              .call(yAxis);       
        
            sites_plotted ++;
        }
    }
    
    svg.attr("height", (dim[1])*sites_plotted + margin.top + margin.bottom);
    
    if (sites_plotted > 0) {
            
        residue_array = [];
        for (k in all_residues) {
            residue_array.push (k);
        }
        all_residues = residue_array;
    
        var legend_dim = {x: width + margin.left, y:20, spacer:15, margin:3, font: 10};
        
        var legend = svg.append("g")
          .attr("class", "legend")
          .attr("x", legend_dim.x)
          .attr("y", legend_dim.y)
          .attr("transform", "translate("+legend_dim.x+","+legend_dim.y+")");
        
        legend.selectAll('g').data(all_residues)
          .enter()
          .append('g')
          .each(function(d, i) {
            var g = d3.select(this);
            g.append("rect")
              .attr("x", legend_dim.spacer)
              .attr("y", i*(legend_dim.spacer + legend_dim.margin))
              .attr("width", legend_dim.spacer)
              .attr("height", legend_dim.spacer)
              .style("fill", _pos_sites_color_map[d]);
    
            g.append("text")
              .attr("x", 2*legend_dim.spacer + legend_dim.font/4)
              .attr("y", (i+1)*(legend_dim.spacer + legend_dim.margin) - legend_dim.margin 
                         - (legend_dim.spacer-legend_dim.font)*2/3)
              .style("fill", _pos_sites_color_map[d])
              .text(function (d) {return d;});
      
          });    
    }

    
}

function _ext_pad (ext) {
    if (ext[1]-ext[0] <= 1) {
        ext[0] -=2;
        ext[1] +=2;
    }
    return ext;
}

function prepare_data_from_keys (data, segments, keys) {

     res = [];
     if (segments.length == 1) {
        for (k in data) {
            if (segments == data[k]["Segment"]) {
                info = [data[k].Date];
                for (k2 in keys) {
                    info.push ((data[k])[keys[k2]]);
                }
                res.push (info);
           }
        }    
    } else {
        info_by_date = {};

        for (k in data) {
            var series_index = segments.indexOf (data[k]["Segment"]);
            if (series_index >= 0) {
                if (!(data[k].Date in info_by_date)) {
                    info_by_date[data[k].Date] = new Array (segments.length + 1);
                     info_by_date[data[k].Date][0] = data[k].Date;
                }
                info_by_date[data[k].Date][series_index + 1] = data[k][keys[0]];
            
           }
        }
        for (k in info_by_date) {
            res.push(info_by_date[k]);
        }   
    }
    
    res.sort (function (a,b) {return a[0]-b[0];});
    return res;
}
 
function plot_pheno_data (data, keys, labels, segment, id, dim) {

    var margin = {top: 20, right: 50, bottom: 60, left: 50},
        width  = dim[0] - margin.left - margin.right,
        height = dim[1] - margin.top - margin.bottom;

    var x = d3.time.scale()
        .range([0, width]);

    var y = d3.scale.linear()
        .range([height, 0]);

    
    var xAxis = d3.svg.axis()
        .scale(x)
        .orient("bottom");

    var yAxis = d3.svg.axis()
        .scale(y)
        .orient("left");

    if (keys.length == 2) {
        var y2 = d3.scale.linear()
            .range([height, 0]);

        var y2Axis = d3.svg.axis()
            .scale(y2)
            .orient("right");
    }
    
    var svg = d3.select("#" + id)
        .attr("width", width + margin.left + margin.right)
        .attr("height", height + margin.top + margin.bottom);

    svg.selectAll("path").remove();
    svg.selectAll("g").remove();    

    sequence_info =  _aux_plot_set_ticks(prepare_data_from_keys (data, segment, keys), x, xAxis);
    
    if (keys.length == 2) {    
        y.domain(_ext_pad (d3.extent(sequence_info, function(d) { return d[1]; }))); 
        y2.domain(_ext_pad (d3.extent(sequence_info, function(d,i) { return d[2]; }))); 
    } else {
        ext = [d3.min (sequence_info.map (function (d) {return d3.min (d.filter (function (e,i) {return i>0;}))})),
               d3.max (sequence_info.map (function (d) {return d3.max (d.filter (function (e,i) {return i>0;}))}))];
               
                
        y.domain(_ext_pad(ext));
    }
    
    
    var plot_svg = svg.append("g")
    .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

      
    var colors = d3.scale.category10();
        
    var line = d3.svg.line()
        .x(function(d) { return x(d.x); })
        .y(function(d) { return y(d.y); });

    plot_svg.append("path")
      .datum(sequence_info.map (function (d) {return {'x' : d[0], 'y': d[1]};}))
      .attr("class", "_evo_line")
      .style('stroke', colors(0))
      .attr("d", line);
      
    
    if (keys.length == 2) {
        var line = d3.svg.line()
            .x(function(d) { return x(d[0]); })
            .y(function(d) { return y2(d[2]); });

        plot_svg.append("path")
          .datum(sequence_info)
          .attr("class", "_evo_line")
          .style('stroke', colors(1))
          .attr("d", line);
    } else {
        for (k = 1; k < segment.length; k++) {
            plot_svg.append("path")
              .datum(sequence_info.map (function (d) {return {'x' : d[0], 'y': d[1+k]};}))
              .attr("class", "_evo_line")
              .style('stroke', colors(k))
              .attr("d", line);
        }
        
        
        for (k = 0; k < segment.length; k++) {
            plot_svg.append ("text")
                .attr ("x", x (sequence_info[sequence_info.length-1][0]) + 10)
                .attr ("y", y (sequence_info[sequence_info.length-1][k+1]))
                .style ("fill", colors (k))
                .text (segment[k])
                .attr ("class", "_evo_legend");
        }
    }


    if (segment.length == 1) {
        var legend_dim = {x: 50, y:20, spacer:25, margin:5, font: 12};
    
        var legend = svg.append("g")
          .attr("class", "_evo_legend")
          .attr("x", legend_dim.x)
          .attr("y", legend_dim.y)
          .attr("transform", "translate("+legend_dim.x+","+legend_dim.y+")");
        
        legend.selectAll('g').data(labels)
          .enter()
          .append('g')
          .each(function(d, i) {
            var g = d3.select(this);
            g.append("rect")
              .attr("x", legend_dim.spacer)
              .attr("y", i*(legend_dim.spacer + legend_dim.margin))
              .attr("width", legend_dim.spacer)
              .attr("height", legend_dim.spacer)
              .style("fill", colors(i));
        
            g.append("text")
              .attr("x", 2*legend_dim.spacer + legend_dim.font/4)
              .attr("y", (i+1)*(legend_dim.spacer + legend_dim.margin) - legend_dim.margin 
                         - (legend_dim.spacer-legend_dim.font)*2/3)
              .style("fill", colors(i))
              .text(function (d) {return d;});
          
          });
    }

    _make_x_axis (plot_svg, xAxis, width, height, labels[0]);
    
    /*var x_g = plot_svg.append("g")
      .attr("class", "x axis")
      .attr("transform", "translate(0," + height + ")")
      .call(xAxis).append("text")
      .attr("transform", "translate("+width+",0)")
      .attr("y", "24")
      .attr("dy", ".71em")
      .style("text-anchor", "end")
      .text("Sample date");*/

    plot_svg.append("g")
      .attr("class", "y axis")
      .call(yAxis)
      .append("text")
      .attr("transform", "rotate(-90)")
      .attr("y", 6)
      .attr("dy", ".71em")
      .style("text-anchor", "end")
      .text(labels[0]);

    if (keys.length == 2) {
        plot_svg.append("g")
          .attr("class", "y axis")
          .attr("transform", "translate("+width + ",0)")
          .call(y2Axis)
          .append("text")
          .attr("transform", "rotate(-90)")
          .attr("y", 6)
          .attr("dy", "-1.5em")
          .attr("dx", -height+10)
          .style("text-anchor", "start")
          .text(labels[1]);
    }

}

function _aux_plot_set_ticks (div_data, x, xAxis) {
    var x_dates = div_data.map (function(d,i) { return d[0]; });
        x_dates.sort();     
        x.domain(d3.extent(x_dates));
        xAxis.tickValues (x_dates).tickFormat (show_date_axis);
    return div_data;
}

function _make_x_axis (plot_svg, xAxis, width, height, label) {
    var x_g = plot_svg.append("g");

    x_g.attr("class", "x axis")
      .attr("transform", "translate(0," + height + ")")
      .call(xAxis);
 
    x_g.selectAll ("text").style("text-anchor", "start").attr ("transform", "rotate(45)").attr("dx","0.5em").attr("dy","0.5em");

    x_g.append("text")
      .attr("transform", "translate("+width+",0)")
      .attr("y", "-24")
      .attr("dy", ".71em")
      .style("text-anchor", "end")
      .text(label);
}

function plot_div_data (data, keys, labels, segment, id, dim, circle_id_in) {
    
    circle_id = [];
    
    if (circle_id_in.length == 1) {
        
        circle_id [1] = keys.indexOf (circle_id_in[0][1]);
        if (circle_id[1] < 0) {
            circle_id = [];
        } else {
            circle_id [0] = circle_id_in[0][0];
        }
    } 
    
    
    
    var margin = {top: 20, right: 50, bottom: 60, left: 50},
        width  = dim[0] - margin.left - margin.right,
        height = dim[1] - margin.top - margin.bottom;
        

    var x = d3.time.scale()
        .range([0, width]);

    var y = d3.scale.linear()
        .range([height, 0]);

    var xAxis = d3.svg.axis()
        .scale(x)
        .orient("bottom");
        //.tickValues (_unique_dates)
        //.tickFormat (function (d) {console.log (d); return d;});

    var yAxis = d3.svg.axis()
        .scale(y)
        .orient("left");

    var svg = d3.select("#" + id)
        .attr("width", width + margin.left + margin.right)
        .attr("height", height + margin.top + margin.bottom);

    svg.selectAll("path").remove();
    svg.selectAll("g").remove();    


    divergence_data = _aux_plot_set_ticks(prepare_data_from_keys (data, segment, keys), x, xAxis);
               
    y.domain([0, d3.max(divergence_data, function(d)  {var values = []; 
        for (k = 1; k < d.length; k++) {
            if (circle_id.length == 2 && circle_id[0] == k  - 1) {
                values.push (d[circle_id[1] + 1] + d[circle_id[0] + 1]/2);
            } else {
                values.push (d[k]);
            }
        }
        return d3.max (values);
    })]);

    
    var plot_svg = svg.append("g")
    .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

      
    var max_x = 0;  
    var circle_gen = d3.svg.symbol().size (function (d) {var r = y(0)-y(d.r); max_x = Math.max (max_x, x(d.x) + r); return r*r;});
      
    var colors = d3.scale.category10();
    
    if (circle_id.length == 2) { 
        
        for (k in divergence_data) {
            plot_svg.append("path")
            .datum({'x' : divergence_data[k][0], 'y': divergence_data[k][circle_id[1]+1], 'r': divergence_data[k][circle_id[0]+1]})
            .attr("transform", function(d) { return "translate(" + x(d.x) + "," + y(d.y) + ")"; })
            .attr("d", circle_gen)
            .attr("class", "_evo_symbol")
            .append("title")
            .text (function (d) { return "Diversity = " + d.r.toFixed (5); });
        }
    
        var overflow = max_x - width;
        if (overflow > 0) {
            svg.attr("width", width + margin.left + margin.right + overflow);
        }
    }
    
    var line = d3.svg.line()
        .x(function(d,i) { return x(d.x); })
        .y(function(d) { return y(d.y); });


    up_to = segment.length == 1 ? keys.length: segment.length;
    
    for (k = 0; k < up_to; k++) {
         if (circle_id.length == 2 && circle_id[0] == k) {
            continue;
         }
         plot_svg.append("path")
          .datum(divergence_data.map (function (d) {return {'x' : d[0], 'y': d[1+k]};}))
          .attr("class", "_evo_line")
          .style('stroke', colors(k))
         .attr("d", line);
   }

    if (segment.length == 1) {
        var legend_dim = {x: 50, y:20, spacer:25, margin:5, font: 12};
    
        var legend = svg.append("g")
          .attr("class", "_evo_legend")
          .attr("x", legend_dim.x)
          .attr("y", legend_dim.y)
          .attr("transform", "translate("+legend_dim.x+","+legend_dim.y+")");
        
        legend.selectAll('g').data(labels)
          .enter()
          .append('g')
          .each(function(d, i0) {
            var g = d3.select(this);
            g.append("rect")
              .attr("x", legend_dim.spacer)
              .attr("y", i0*(legend_dim.spacer + legend_dim.margin))
              .attr("width", legend_dim.spacer)
              .attr("height", legend_dim.spacer)
              .style("fill", function () {if (i0 == circle_id[1] + 1) return "#DDDDDD"; return colors(i0);});
        
            g.append("text")
              .attr("x", 2*legend_dim.spacer + legend_dim.font/4)
              .attr("y", (i0+1)*(legend_dim.spacer + legend_dim.margin) - legend_dim.margin 
                         - (legend_dim.spacer-legend_dim.font)*2/3)
              .style("fill", function () {if (i0 == circle_id[1] + 1) return "#DDDDDD"; return colors(i0);})
              .text(function (d) {return d;});
          
          });
    } else {
        for (k = 0; k < segment.length; k++) {
                plot_svg.append ("text")
                    .attr ("x", x (divergence_data[divergence_data.length-1][0]) + 10)
                    .attr ("y", y (divergence_data[divergence_data.length-1][k+1]))
                    .style ("fill", colors (k))
                    .text (segment[k])
                    .attr ("class", "_evo_legend");
        }
    }


    _make_x_axis (plot_svg, xAxis, width, height, "Sample Date");
    
    

    plot_svg.append("g")
      .attr("class", "y axis")
      .call(yAxis)
      .append("text")
      .attr("transform", "rotate(-90)")
      .attr("y", 6)
      .attr("dy", ".71em")
      .style("text-anchor", "end")
      .text(segment.length == 1 ? "Expected substitutions per site" : labels[0]);
}

//helpers

function _add_dict_values (d, excl) {
    var sum = 0.;
    for (var key in d) {
        if (excl && key in excl) {
            continue;
        }
        sum += d[key];
    }
    return sum;
}

function draw_positional_data () {
    if ($( "#_pos_entropy" ).prop( "checked" )) {
        handle_positional_data ([_pos_overall_data[0], _pos_overall_data[1].map (function (d) { return [d[4]];})], _pos_data.map (function (d) { return [d[0], d[1].map (function (d) {return [d[4]]})];}),  "_pos_main_plot", [850,_positional_data_plot_height], ["Entropy"]);
    } else {
        handle_positional_data ([_pos_overall_data[0], _pos_overall_data[1].map (function (d) { return [d[0],d[1]];})], _pos_data.map (function (d) { return [d[0], d[1].map (function (d) {return [d[0],d[1]]})];}),  "_pos_main_plot", [850,_positional_data_plot_height], ["Mean dS", "Mean dN"]);
    }
}

function check_positive_selection (mx) {
    return mx.map (function (d, i) {return [i, d[2]]; }).filter (function (d) {return d [1] >= 0.95;}).map (function (d) {return d[0] + 1;});    
}

function get_site_residues (data, site) {
    var all_residues = {};
    
    for (k in data[site]) {
        if (k != "HXB2") {
            for (r in data[site][k]) {
                all_residues[r] = 1;
            }
        }
    }
    
    return d3.keys (all_residues).sort();
}      

function  get_db_residue_counts (index, subtype) {
    if (_db_frequency_information) {
        if (subtype == undefined || subtype == 'Combined') {
            pull_from  = _db_frequency_information["overall"];
        } else {
            pull_from = _db_frequency_information["subtype"][subtype];
        }
        if (index in pull_from) {
            var residues = pull_from[index];
            return residues;
        }
    }
    
    return null;
}

function add_options_to_select (id, options) {
    d3.select ("#" + id)
              .selectAll ("option")
              .data (options)
              .enter()
              .append ("option")
              .text (function (d) {return d[0];})
              .attr ('value', function (d) {return d[1];});
              
}