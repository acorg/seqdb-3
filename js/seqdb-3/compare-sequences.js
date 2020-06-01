
function main()
{
    console.log("compare_sequences_data", compare_sequences_data);
    show_most_frequent_per_group(document.querySelector("#compare-sequences .most-frequent-per-group"))
    show_frequency_per_group(document.querySelector("#compare-sequences .frequency-per-group"))
    show_full_sequences(document.querySelector("#compare-sequences .full-sequences"));
}

// --------------------------------------------------------------------------------

function show_full_sequences(div)
{
    const title = document.createElement("p");
    title.classList.add("title");
    title.innerHTML = "Full sequences";
    div.appendChild(title);
    
    const add_seqence = function(tr, seq_s) {
        const seq = [...seq_s];
        seq.forEach(function(aa, pos0) {
            const aa_td = document.createElement("td");
            aa_td.classList.add(`aa${aa}`);
            aa_td.classList.add("aa");
            if (pos0 % 10 == 9)
                aa_td.classList.add("sep-left-zero");
            else if (pos0 % 10 == 4)
                aa_td.classList.add("sep-left-six");
            aa_td.innerHTML = aa;
            tr.appendChild(aa_td);
        });
    };

    const add_ruler = function() {
        const tr = document.createElement("tr");
        tr.classList.add("aa-ruler");
        tr.innerHTML = `<td colspan="2"></td>`;
        for (let pos1 = 1; pos1 < 550; ++pos1) {
            const td = document.createElement("td");
            if (pos1 % 10 == 0) {
                td.classList.add("sep-left-zero");
                let pos_s = "" + pos1;
                td.innerHTML = pos_s;
                td.setAttribute("colspan", pos_s.length);
                pos1 += pos_s.length - 1;
            }
            else if (pos1 % 10 == 5) {
                td.classList.add("sep-left-six");
                td.innerHTML = ".";
            }
            else
                td.innerHTML = ".";
            tr.appendChild(td);
        }
        return tr;
    };

    const tab1 = document.createElement("table");
    for (let group of compare_sequences_data.groups) {
        tab1.appendChild(add_ruler());
        group.seq.forEach(function(id_seq, index) {
            const tr = document.createElement("tr");
            if (index == 0)
                tr.innerHTML = `<td class="group-name" rowspan="${group.seq.length}">${group.name}</td><td class="seq-id">${id_seq.id}</td>`;
            else
                tr.innerHTML = `<td class="seq-id">${id_seq.id}</td>`;
            add_seqence(tr, id_seq.seq);
            tab1.appendChild(tr);
        });
    }
    tab1.appendChild(add_ruler());
    div.appendChild(tab1);
}

// --------------------------------------------------------------------------------

function show_most_frequent_per_group(div)
{
    const title = document.createElement("p");
    title.classList.add("title");
    title.innerHTML = "Most frequent per group";
    div.appendChild(title);
}

// --------------------------------------------------------------------------------

function show_frequency_per_group(div)
{
    const title = document.createElement("p");
    title.classList.add("title");
    title.innerHTML = "Frequency per group";
    div.appendChild(title);
}

// --------------------------------------------------------------------------------

document.addEventListener("DOMContentLoaded", main);
